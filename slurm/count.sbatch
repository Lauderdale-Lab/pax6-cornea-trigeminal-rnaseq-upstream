#!/usr/bin/env bash
#
# slurm/count.sbatch -- one gene-count matrix per study, in one featureCounts
# invocation over every library of every dataset in that study.
#
# Refuses to count unless, for every dataset in the study:
#   * every library in samples.tsv has a verified BAM,
#   * the strand check in THIS run agrees with config/datasets.tsv,
#   * all datasets share one strand setting and one layout.
#
#SBATCH --job-name=pax6_count
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00

set -euo pipefail
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
# shellcheck source=../lib/provenance.sh
source "$PAX6_HOME/lib/provenance.sh"
load_run "$PAX6_RUN"
STUDY="$PAX6_STUDY"
require_study "$STUDY"
OUT=$(counts_dir "$STUDY")
[[ -e "$OUT/featureCounts.txt" ]] && die "$OUT/featureCounts.txt exists. Counts are never overwritten; start a new run."
mkdir -p "$OUT"
load_modules Subread SAMtools
record_versions count Subread SAMtools

# --- every dataset: measured == recorded, and all agree ---------------------
STRAND=""; LAYOUT=""; BAMS=(); : > "$OUT/sample_map.tsv.tmp"
while read -r DS; do
  R="$(strand_dir "$DS")/RESULT.tsv"
  [[ -s "$R" ]] || die "$DS: no strand result in this run ($R). Run: bin/pax6 submit $RUN $DS"
  IFS=$'\t' read -r ML MS _ < <(sed -n 2p "$R")
  RL=$(dataset_field "$DS" layout); RS=$(dataset_field "$DS" strand)
  [[ "$ML" == "$RL" && "$MS" == "$RS" ]] || die "$DS: measured layout=$ML strand=$MS but config/datasets.tsv says layout=$RL strand=$RS.
        Correct config/datasets.tsv (and commit it), or investigate. Do not count until they agree."
  [[ -z "$STRAND" || "$STRAND" == "$MS" ]] || die "Datasets of '$STUDY' differ in strandedness ($STRAND vs $MS). One featureCounts run applies one -s; count them as separate studies."
  [[ -z "$LAYOUT" || "$LAYOUT" == "$ML" ]] || die "Datasets of '$STUDY' mix PE and SE. Fragment and read counts are not on one scale; split the study."
  STRAND="$MS"; LAYOUT="$ML"

  BATCH=$(dataset_field "$DS" batch); ALN=$(align_dir "$DS")
  while IFS=$'\t' read -r S _ _; do
    B="$ALN/$S.bam"
    bam_ok "$B" || die "$DS/$S: BAM missing, unindexed or truncated: $B"
    BAMS+=("$B")
    printf '%s\t%s\t%s\t%s\t%s\n' "$S" "$STUDY" "$DS" "$BATCH" "$B" >> "$OUT/sample_map.tsv.tmp"
  done < "$(samples_tsv "$DS")"
  log "$DS: $(wc -l < "$(samples_tsv "$DS")") libraries, layout $ML, strand -s $MS"
done < <(study_datasets "$STUDY")

{ printf 'sample_id\tstudy\tdataset\tbatch\tbam\n'; cat "$OUT/sample_map.tsv.tmp"; } > "$OUT/sample_map.tsv"
rm -f "$OUT/sample_map.tsv.tmp"
DUPS=$(tail -n +2 "$OUT/sample_map.tsv" | cut -f1 | sort | uniq -d)
[[ -z "$DUPS" ]] || die "Sample IDs repeated within study '$STUDY': $DUPS. Two columns mapping to one metadata row is silent corruption."
[[ -n "${PAX6_EXPECT_N:-}" && "${#BAMS[@]}" != "$PAX6_EXPECT_N" ]] && die "Expected $PAX6_EXPECT_N libraries, found ${#BAMS[@]}."

# --- count ---------------------------------------------------------------
PE=(); [[ "$LAYOUT" == PE ]] && read -r -a PE <<< "$FC_PE_ARGS"
# shellcheck disable=SC2206
FC=($FC_ARGS)
CMD=(featureCounts -T "$THREADS" -a "$ANNOTATION" -F GTF "${PE[@]}" "${FC[@]}" -s "$STRAND" -o "$OUT/featureCounts.txt" "${BAMS[@]}")
log "counting ${#BAMS[@]} libraries"
"${CMD[@]}" 2> "$OUT/featureCounts.log"

# --- tidy matrix: gene_id x sample_id, no comment line, no coordinates ----
awk -F'\t' 'NR == 1 { next }
  NR == 2 { printf "gene_id"; for (i = 7; i <= NF; i++) { n = split($i, p, "/"); s = p[n]; sub(/\.bam$/, "", s); printf "\t%s", s } ; print ""; next }
  { printf "%s", $1; for (i = 7; i <= NF; i++) printf "\t%s", $i; print "" }' \
  "$OUT/featureCounts.txt" > "$OUT/gene_counts.tsv"

# --- provenance --------------------------------------------------------------
{
  printf 'key\tvalue\n'
  printf 'run\t%s\n' "$RUN"
  printf 'study\t%s\n' "$STUDY"
  printf 'created\t%s\n' "$(date -Iseconds)"
  printf 'pipeline_commit\t%s\n' "$(run_info "$RUN" pipeline_commit)"
  printf 'assembly\t%s\n' "$ASSEMBLY"
  printf 'annotation\t%s\n' "$ANNOTATION"
  printf 'index\t%s (%s)\n' "$HISAT2_INDEX" "$INDEX_TYPE"
  printf 'layout\t%s\n' "$LAYOUT"
  printf 'strand\t%s\n' "$STRAND"
  printf 'libraries\t%s\n' "${#BAMS[@]}"
  printf 'featurecounts\t%s\n' "$(featureCounts -v 2>&1 | grep -o 'v[0-9.]*' | head -n1)"
  printf 'command\t%s\n' "${CMD[*]:0:$((${#CMD[@]} - ${#BAMS[@]}))} <${#BAMS[@]} BAMs>"
  printf 'matrix_md5\t%s\n' "$(md5sum < "$OUT/gene_counts.tsv" | cut -d' ' -f1)"
} > "$OUT/PROVENANCE.tsv"

# --- quality gate ------------------------------------------------------------
assignment_table "$OUT/featureCounts.txt.summary" > "$OUT/assignment_qc.tsv"
print_table "$OUT/assignment_qc.tsv"
write_provenance "$STUDY" "$OUT"
log "provenance: $OUT/provenance/   (Methods draft: provenance/METHODS.md)"
# Runs are processed on scratch, which is purged. Keep the matrix, its
# provenance and the run's QC and records on /work now, before the quality
# gate below, so even a flagged matrix is preserved for inspection.
KEPT=$(keep_results "$RUN")
log "kept on /work: $KEPT/counts/$STUDY"
LOW=$(awk -F'\t' -v m="$MIN_ASSIGN_RATIO" 'NR > 1 && $4 < m { print $1 }' "$OUT/assignment_qc.tsv")
log "matrix: $OUT/gene_counts.tsv"
if [[ -n "$LOW" ]]; then
  warn "Assigned/(Assigned+NoFeatures) below ${MIN_ASSIGN_RATIO}% for: $(tr '\n' ' ' <<<"$LOW")"
  warn "Near 50% across the board means the strand setting is wrong; one low library is more likely degraded or low-complexity RNA. Resolve before using this matrix."
  exit 3
fi
