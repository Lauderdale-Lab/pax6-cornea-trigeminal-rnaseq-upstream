#!/usr/bin/env bash
#
# align_batch.sh -- align one batch to sorted, indexed BAMs.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# USAGE
#   sbatch align_batch.sh <study> <batch>
#   sbatch --array=1-21%6 align_batch.sh <study> <batch>
#
#   Assembly defaults to the first row of assemblies.tsv. To use another:
#     sbatch --export=ALL,PAX6_ASSEMBLY=GRCm39 align_batch.sh <study> <batch>
#
#   BAMs are written to <mapping_root>/<assembly>/, so the same libraries can
#   be aligned to several assemblies without colliding, and a count matrix can
#   never mix them.
#
# SCOPE
#   Produces BAMs and nothing else. StringTie and prepDE were removed: their
#   output was never consumed by the analysis, both were mis-parameterised
#   (--rf on unstranded libraries; prepDE at its default 75 bp against 150 bp
#   reads), and wrong files left in an output directory invite a later script
#   to pick them up.
#
#   Counting is NOT here. Counting spans batches by definition, and keeping it
#   inside a per-batch script is what produced two batches counted months
#   apart in separate jobs. See count_study.sh.
#
#SBATCH --job-name=align_batch
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=align_%A_%a.out
#SBATCH --error=align_%A_%a.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$PIPELINE_DIR/pipeline_common.sh"

STUDY="${1:-}"; BATCH="${2:-}"
[[ -n "$STUDY" && -n "$BATCH" ]] || die "Usage: sbatch align_batch.sh <study> <batch>"

require_manifest; require_study "$STUDY"
TRIMMED_DIR=$(manifest_field "$STUDY" "$BATCH" 5)
MAPPING_DIR=$(mapping_dir_for "$STUDY" "$BATCH")
[[ -n "$TRIMMED_DIR" ]] || die "No manifest row for $STUDY/$BATCH"
[[ -d "$TRIMMED_DIR" ]] || die "Trimmed FASTQ directory not found: $TRIMMED_DIR"
[[ -f "$ANNOTATION" ]]  || die "Annotation not found: $ANNOTATION"
[[ -s "${HISAT2_INDEX}.1.ht2" ]] || die "HISAT2 index not found: ${HISAT2_INDEX}.1.ht2"

load_pinned HISAT2 SAMtools
log "$STUDY/$BATCH   assembly $ASSEMBLY ($INDEX_TYPE index)"
log "  annotation : $ANNOTATION"
log "  index      : $HISAT2_INDEX"
log "  BAMs       : $MAPPING_DIR"
report_versions
mkdir -p "$MAPPING_DIR"

SPLICE="$MAPPING_DIR/splice_sites.txt"
if [[ ! -s "$SPLICE" ]]; then
  SPLICE_PY=$(find "$(dirname "$(command -v hisat2)")/.." -name hisat2_extract_splice_sites.py 2>/dev/null | head -n1)
  [[ -n "$SPLICE_PY" ]] || die "hisat2_extract_splice_sites.py not found beside hisat2"
  log "extracting splice sites"
  python3 "$SPLICE_PY" "$ANNOTATION" > "$SPLICE"
  [[ -s "$SPLICE" ]] || die "Splice-site extraction produced an empty file"
fi

# Layout from the files present, not from the manifest.
mapfile -t R1S < <(find_r1 "$TRIMMED_DIR")
mapfile -t SES < <(find_se "$TRIMMED_DIR")
if   (( ${#R1S[@]} > 0 )); then LAYOUT=PE; INPUTS=("${R1S[@]}")
elif (( ${#SES[@]} > 0 )); then LAYOUT=SE; INPUTS=("${SES[@]}")
else die "No trimmed FASTQs in $TRIMMED_DIR"; fi
log "$LAYOUT, ${#INPUTS[@]} sample(s)"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  (( idx >= 0 && idx < ${#INPUTS[@]} )) || die "Array index out of range (1-${#INPUTS[@]})"
  INPUTS=("${INPUTS[$idx]}")
fi

n_done=0; n_skip=0
for IN in "${INPUTS[@]}"; do
  S=$(sample_from_r1 "$IN")
  BAM="$MAPPING_DIR/${S}_hisat2.sorted.bam"
  SUMMARY="$MAPPING_DIR/${S}_hisat2_summary.txt"

  # Resume on a VERIFIED product, not on a filename. The old pipeline skipped
  # whenever a path existed, so a BAM truncated by a killed job read as done.
  if bam_ok "$BAM"; then log "$S: complete BAM present, skipping"; n_skip=$((n_skip+1)); continue; fi
  [[ -e "$BAM" ]] && { warn "$S: incomplete BAM, redoing"; rm -f "$BAM" "$BAM.bai"; }

  log "$S: aligning ($LAYOUT)"
  # Piped into samtools sort: the intermediate SAM is never written. At this
  # depth each one runs to roughly 200 GB, and the old pipeline kept them all.
  if [[ "$LAYOUT" == PE ]]; then
    R2=$(mate_of "$IN") || die "No mate for $IN"
    hisat2 -p "$THREADS" --dta --known-splicesite-infile "$SPLICE" \
        --summary-file "$SUMMARY" -x "$HISAT2_INDEX" -1 "$IN" -2 "$R2" \
      | samtools sort -@ 4 -m 2G -o "$BAM" -
  else
    hisat2 -p "$THREADS" --dta --known-splicesite-infile "$SPLICE" \
        --summary-file "$SUMMARY" -x "$HISAT2_INDEX" -U "$IN" \
      | samtools sort -@ 4 -m 2G -o "$BAM" -
  fi

  samtools index -@ 4 "$BAM"
  bam_ok "$BAM" || die "$S: BAM failed verification after alignment"
  n_done=$((n_done+1))
done

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  REPORT="$MAPPING_DIR/AlignmentSummary_${STUDY}_${BATCH}_${ASSEMBLY}.tsv"
  {
    printf 'sample_id\ttotal\toverall_alignment_rate\tconcordant_exactly_once_pct\n'
    for F in "$MAPPING_DIR"/*_hisat2_summary.txt; do
      [[ -e "$F" ]] || continue
      awk -v s="$(basename "${F%_hisat2_summary.txt}")" '
        /reads; of these:/ { total = $1 }
        /aligned concordantly exactly 1 time/ { gsub(/[()%]/, "", $2); conc = $2 }
        /overall alignment rate/ { gsub(/%/, "", $1); rate = $1 }
        END { printf "%s\t%s\t%s\t%s\n", s, total, rate, (conc == "" ? "NA" : conc) }' "$F"
    done
  } > "$REPORT"
  log "alignment summary: $REPORT"
fi

log "$STUDY/$BATCH: $n_done aligned, $n_skip already complete"
log "NEXT: sbatch strand_check.sh $STUDY $BATCH   (assembly $ASSEMBLY)"
