#!/usr/bin/env bash
#
# slurm/align.sbatch -- one library: HISAT2 -> sorted, indexed BAM.
# Array task N handles line N of the dataset's samples.tsv, reading the
# trimmed reads that qc_trim.sbatch wrote.
#
#SBATCH --job-name=pax6_align
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
load_run "$PAX6_RUN"
DATASET="$PAX6_DATASET"
IFS=$'\t' read -r SAMPLE _R1 R2 < <(sample_line "$(samples_tsv "$DATASET")" "${SLURM_ARRAY_TASK_ID:?not an array task}")

TRIM=$(trim_dir "$DATASET"); ALN=$(align_dir "$DATASET")
mkdir -p "$ALN"
BAM="$ALN/${SAMPLE}.bam"
log "$RUN / $DATASET / $SAMPLE  ->  $ASSEMBLY ($INDEX_TYPE index)"

load_modules HISAT2 SAMtools
if bam_ok "$BAM"; then log "BAM present and verified; nothing to do"; exit 0; fi
rm -f "$BAM" "$BAM.bai"

[[ -s "${HISAT2_INDEX}.1.ht2" ]] || die "HISAT2 index not found: ${HISAT2_INDEX}.1.ht2"
[[ -s "$ANNOTATION" ]]           || die "Annotation not found: $ANNOTATION"
[[ -f "$TRIM/$SAMPLE.done" ]]    || die "Trimming not complete for $SAMPLE (no $TRIM/$SAMPLE.done)"

# Known splice sites from the annotation, extracted once per run. Array tasks
# start together, so the first to take the lock writes the file.
SPLICE="$RUN_DIR/reference/splice_sites.txt"
mkdir -p "$RUN_DIR/reference"
(
  flock 9
  if [[ ! -s "$SPLICE" ]]; then
    PY=$(find "$(dirname "$(command -v hisat2)")/.." -name hisat2_extract_splice_sites.py 2>/dev/null | head -n1)
    [[ -n "$PY" ]] || die "hisat2_extract_splice_sites.py not found beside hisat2"
    python3 "$PY" "$ANNOTATION" > "$SPLICE.tmp" && mv "$SPLICE.tmp" "$SPLICE"
    log "splice sites: $(wc -l < "$SPLICE") from $(basename "$ANNOTATION")"
  fi
) 9> "$SPLICE.lock"
[[ -s "$SPLICE" ]] || die "Splice-site file is empty: $SPLICE"

if [[ "$R2" != - ]]; then READS=(-1 "$TRIM/${SAMPLE}_1.fq.gz" -2 "$TRIM/${SAMPLE}_2.fq.gz")
else                      READS=(-U "$TRIM/${SAMPLE}.fq.gz"); fi

# hisat2 streams into samtools sort, so no SAM is ever written; the BAM gets
# its final name only after it verifies.
# shellcheck disable=SC2086  # HISAT2_ARGS is a deliberate word list
hisat2 -p "$THREADS" $HISAT2_ARGS --known-splicesite-infile "$SPLICE" \
    --summary-file "$ALN/${SAMPLE}.hisat2.txt" --new-summary \
    -x "$HISAT2_INDEX" "${READS[@]}" \
  | samtools sort -@ 4 -m 1500M -O bam -T "${TMPDIR:-$ALN}/${SAMPLE}.sorttmp" -o "$BAM.tmp" -
mv "$BAM.tmp" "$BAM"
samtools index -@ 4 "$BAM"
bam_ok "$BAM" || die "BAM failed verification after alignment: $BAM"

grep -iE 'overall alignment rate|Overall alignment rate' "$ALN/${SAMPLE}.hisat2.txt" || true
log "done"
