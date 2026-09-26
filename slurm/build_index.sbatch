#!/usr/bin/env bash
#
# slurm/build_index.sbatch -- build a HISAT2 index for a genome build listed
# in config/assemblies.tsv. Submitted by `bin/pax6 index <assembly>`.
#
# The genome FASTA and matching annotation must already be in
# $PAX6_ROOT/reference/<assembly>/, uncompressed.
#
# index_type plain   : ~8 GB RAM, a few hours. Splice awareness still comes
#                      from --known-splicesite-infile at alignment time.
# index_type ss_exon : splice sites and exons built in; needs ~200 GB RAM
#                      (submit with --partition=highmem_p --mem=200G).
# snp_tran is not built here. If you ever build one, confirm afterwards with
#   hisat2-inspect --snp <prefix> | wc -l
# The GRCm38 "genome_snp_tran" index turned out to contain 0 SNPs because its
# VCF named chromosomes 1, 2, ... and its FASTA chr1, chr2, ...
#
#SBATCH --job-name=pax6_index
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00

set -euo pipefail
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
resolve_assembly "$PAX6_ASSEMBLY"
load_modules HISAT2

[[ -s "$GENOME_FASTA" ]] || die "Genome FASTA not found or empty: $GENOME_FASTA"
[[ -s "$ANNOTATION" ]]   || die "Annotation not found or empty: $ANNOTATION"
[[ -s "${HISAT2_INDEX}.1.ht2" ]] && { log "Index already present: ${HISAT2_INDEX}. Delete its .ht2 files to rebuild."; exit 0; }
mkdir -p "$(dirname "$HISAT2_INDEX")"

# Annotation and genome must use the same sequence names. A mismatch builds
# without error and later aligns fine while counting assigns almost nothing.
SHARED=$(comm -12 <(grep -v '^#' "$ANNOTATION" | cut -f1 | sort -u) \
                  <(grep '^>' "$GENOME_FASTA" | sed 's/^>//; s/[[:space:]].*//' | sort -u) | wc -l)
(( SHARED >= 10 )) || die "Only $SHARED sequence names shared by $(basename "$ANNOTATION") and $(basename "$GENOME_FASTA"). Different builds, or chr1 vs 1 naming."
log "$SHARED sequence names shared by genome and annotation"

HELPER() { find "$(dirname "$(command -v hisat2)")/.." -name "$1" | head -n1; }
case "$INDEX_TYPE" in
  plain)
    hisat2-build -p "$THREADS" "$GENOME_FASTA" "$HISAT2_INDEX" ;;
  ss_exon)
    SS="$GENOME_DIR/splice_sites.txt"; EX="$GENOME_DIR/exons.txt"
    python3 "$(HELPER hisat2_extract_splice_sites.py)" "$ANNOTATION" > "$SS"
    python3 "$(HELPER hisat2_extract_exons.py)"        "$ANNOTATION" > "$EX"
    hisat2-build -p "$THREADS" --ss "$SS" --exon "$EX" "$GENOME_FASTA" "$HISAT2_INDEX" ;;
  *)
    die "index_type '$INDEX_TYPE' is not built by this script (use plain or ss_exon)." ;;
esac

for i in 1 2 3 4 5 6 7 8; do [[ -s "${HISAT2_INDEX}.$i.ht2" ]] || die "Build finished but ${HISAT2_INDEX}.$i.ht2 is missing"; done
log "index complete: $HISAT2_INDEX"
log "splice sites built in: $(hisat2-inspect --ss "$HISAT2_INDEX" 2>/dev/null | wc -l)   SNPs: $(hisat2-inspect --snp "$HISAT2_INDEX" 2>/dev/null | wc -l)"
