#!/usr/bin/env bash
#
# build_hisat2_index.sh -- build a HISAT2 index for an assembly.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# USAGE
#   sbatch build_hisat2_index.sh <assembly>          # e.g. GRCm39
#
# The assembly must have a row in assemblies.tsv, and its genome FASTA and
# annotation must already be in that row's ref_dir, uncompressed.
#
# READ THIS BEFORE BUILDING
#
#   You probably do not need a transcript-aware index.
#     The `--ss/--exon` build bakes splice sites into the index and needs
#     roughly 200 GB of RAM and many hours. But align_batch.sh already passes
#     `--known-splicesite-infile`, which supplies the same splice sites AT
#     ALIGNMENT TIME for no build cost at all. A plain index (~8 GB RAM, a few
#     hours) plus that flag gives splice-aware alignment.
#
#   The GRCm38 index is transcript-aware but NOT SNP-aware.
#     Its directory is named `genome_snp_tran`, but hisat2-inspect finds 0 SNPs
#     and 274,187 splice sites in it. The VCF used at build time named
#     chromosomes 1, 2, ... while the FASTA used chr1, chr2, ..., so every
#     variant was skipped without an error. If you build a SNP-aware index,
#     confirm afterwards with `hisat2-inspect --snp <prefix> | wc -l`.
#
#   Prebuilt indexes use Ensembl chromosome names.
#     The HISAT2 authors distribute prebuilt indexes, but theirs use Ensembl
#     naming (1, 2, ...), which does not match GENCODE annotation (chr1, chr2,
#     ...). A mismatch lets alignment succeed while counting assigns almost
#     nothing. Check the names in any index against the annotation before
#     using it.
#
#   Set INDEX_TYPE in assemblies.tsv to whatever you actually build, and put
#     it in the Methods. "Aligned to GRCm39" is not a sufficient description.
#
#SBATCH --job-name=hisat2_build
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=hisat2_build_%j.out
#SBATCH --error=hisat2_build_%j.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"

TARGET="${1:-}"
[[ -n "$TARGET" ]] || { echo "Usage: sbatch build_hisat2_index.sh <assembly>" >&2; exit 1; }
export PAX6_ASSEMBLY="$TARGET"
source "$PIPELINE_DIR/pipeline_common.sh"

load_pinned HISAT2
report_versions

log "assembly    : $ASSEMBLY"
log "genome      : $GENOME_FASTA"
log "annotation  : $ANNOTATION"
log "index       : $HISAT2_INDEX  (declared type: $INDEX_TYPE)"

[[ -s "$GENOME_FASTA" ]] || die "Genome FASTA not found or empty: $GENOME_FASTA
        Download the primary assembly FASTA for $ASSEMBLY into $GENOME_DIR and gunzip it."
[[ -s "$ANNOTATION" ]]   || die "Annotation not found or empty: $ANNOTATION
        It MUST be the GENCODE release built on $ASSEMBLY. vM25 is the last
        release on GRCm38; anything later is on GRCm39."

if [[ -s "${HISAT2_INDEX}.1.ht2" ]]; then
  log "Index already present: ${HISAT2_INDEX}.1.ht2 -- nothing to do."
  log "Delete the .ht2 files first if you intend to rebuild."
  exit 0
fi
mkdir -p "$(dirname "$HISAT2_INDEX")"

# Sanity check that annotation and genome refer to the same sequences. A vM25
# GTF against a GRCm39 FASTA parses without error and produces coordinates that
# do not correspond to genes -- silent, and catastrophic.
log "checking that the annotation matches the genome"
GTF_CHR=$(awk '$0 !~ /^#/ { print $1 }' "$ANNOTATION" | sort -u | head -30)
FA_CHR=$(grep '^>' "$GENOME_FASTA" | sed 's/^>//; s/ .*//' | sort -u | head -60)
OVERLAP=$(comm -12 <(printf '%s\n' "$GTF_CHR" | sort) <(printf '%s\n' "$FA_CHR" | sort) | wc -l)
(( OVERLAP >= 10 )) || die "Only $OVERLAP sequence names shared between the annotation and the genome.
        These files are probably from different assemblies, or use different
        naming (chr1 vs 1). Resolve this before building."
log "  $OVERLAP shared sequence names -- consistent"

case "$INDEX_TYPE" in
  tran|plain)
    # Plain index. Splice awareness comes from --known-splicesite-infile at
    # alignment time, which align_batch.sh always passes.
    log "building plain index (splice sites supplied at alignment time)"
    hisat2-build -p "$THREADS" "$GENOME_FASTA" "$HISAT2_INDEX"
    ;;
  ss_exon)
    # Transcript-aware index. Needs a high-memory node; check that the
    # partition and --mem above are adequate before submitting.
    SS="$GENOME_DIR/splice_sites.txt"; EX="$GENOME_DIR/exons.txt"
    SS_PY=$(find "$(dirname "$(command -v hisat2)")/.." -name hisat2_extract_splice_sites.py | head -n1)
    EX_PY=$(find "$(dirname "$(command -v hisat2)")/.." -name hisat2_extract_exons.py | head -n1)
    [[ -n "$SS_PY" && -n "$EX_PY" ]] || die "HISAT2 helper scripts not found"
    [[ -s "$SS" ]] || python3 "$SS_PY" "$ANNOTATION" > "$SS"
    [[ -s "$EX" ]] || python3 "$EX_PY" "$ANNOTATION" > "$EX"
    log "building transcript-aware index -- this needs ~200 GB RAM"
    hisat2-build -p "$THREADS" --ss "$SS" --exon "$EX" "$GENOME_FASTA" "$HISAT2_INDEX"
    ;;
  snp_tran)
    die "snp_tran indexes are not built here. They require a variant VCF whose
        chromosome names match the FASTA, and substantial memory.
        Download a prebuilt snp_tran index, or set index_type to 'tran' in
        assemblies.tsv and accept a plain index -- recording the difference."
    ;;
  *)
    die "Unknown index_type '$INDEX_TYPE' in assemblies.tsv. Use tran or ss_exon."
    ;;
esac

MISSING=0
for i in {1..8}; do [[ -s "${HISAT2_INDEX}.${i}.ht2" ]] || MISSING=1; done
(( MISSING == 0 )) || die "Index build finished but not all 8 .ht2 files are present."

log "index complete: $HISAT2_INDEX"
log "NEXT: PAX6_ASSEMBLY=$ASSEMBLY sbatch align_batch.sh <study> <batch>"
log "Record in Methods: assembly $ASSEMBLY, annotation $(basename "$ANNOTATION"), index type $INDEX_TYPE."
