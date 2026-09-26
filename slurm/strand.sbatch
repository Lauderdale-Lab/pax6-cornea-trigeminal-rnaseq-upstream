#!/usr/bin/env bash
#
# slurm/strand.sbatch -- measure layout and strandedness for one dataset.
# featureCounts at -s 0, 1 and 2 on a few BAMs. Writes strand/RESULT.tsv,
# which count.sbatch compares with config/datasets.tsv before counting.
#
# Why this runs every time: both Lauderdale batches were once counted at -s 2
# on the assumption of a dUTP library. They are unstranded, and half of every
# library's reads were discarded without any error.
#
#SBATCH --job-name=pax6_strand
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
load_run "$PAX6_RUN"
DATASET="$PAX6_DATASET"
ALN=$(align_dir "$DATASET"); OUT=$(strand_dir "$DATASET")
mkdir -p "$OUT"
load_modules Subread SAMtools
record_versions strand Subread SAMtools

mapfile -t BAMS < <(find "$ALN" -maxdepth 1 -name '*.bam' | sort | head -n "$STRAND_TEST_SAMPLES")
(( ${#BAMS[@]} > 0 )) || die "No BAMs in $ALN"

LAYOUT=SE
[[ "$(samtools view -c -f 1 "${BAMS[0]}")" -gt 0 ]] && LAYOUT=PE
PE=(); [[ "$LAYOUT" == PE ]] && read -r -a PE <<< "$FC_PE_ARGS"
# shellcheck disable=SC2206
FC=($FC_ARGS)

TABLE="$OUT/strand_table.tsv"
printf 'sample_id\ts\tassigned_fraction\n' > "$TABLE"
for BAM in "${BAMS[@]}"; do
  bam_ok "$BAM" || die "BAM failed verification: $BAM"
  S=$(basename "$BAM" .bam)
  for X in 0 1 2; do
    featureCounts -T "$THREADS" -a "$ANNOTATION" -F GTF "${PE[@]}" "${FC[@]}" \
      -s "$X" -o "$OUT/${S}_s${X}.txt" "$BAM" > "$OUT/${S}_s${X}.log" 2>&1
    printf '%s\t%s\t%s\n' "$S" "$X" "$(assigned_fraction "$OUT/${S}_s${X}.txt.summary")" >> "$TABLE"
  done
done

read -r A0 A1 A2 < <(awk -F'\t' 'NR > 1 { s[$2] += $3; n[$2]++ } END { printf "%f %f %f\n", s[0]/n[0], s[1]/n[1], s[2]/n[2] }' "$TABLE")
CALL=$(strand_call "$A0" "$A1" "$A2")
printf 'layout\tstrand\tassigned_s0\tassigned_s1\tassigned_s2\tlibraries_tested\n%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$LAYOUT" "$CALL" "$A0" "$A1" "$A2" "${#BAMS[@]}" > "$OUT/RESULT.tsv"

print_table "$TABLE"
EXPECT_L=$(dataset_field "$DATASET" layout); EXPECT_S=$(dataset_field "$DATASET" strand)
echo
echo "MEASURED   layout=$LAYOUT strand=$CALL"
echo "RECORDED   layout=$EXPECT_L strand=$EXPECT_S   (config/datasets.tsv)"
if [[ "$CALL" == AMBIGUOUS ]]; then
  echo "The three settings match no expected pattern. Inspect $OUT before going further."
elif [[ "$LAYOUT" != "$EXPECT_L" || "$CALL" != "$EXPECT_S" ]]; then
  echo "DISAGREE: record the measured values in config/datasets.tsv (and commit) before counting."
else
  echo "AGREE"
fi
