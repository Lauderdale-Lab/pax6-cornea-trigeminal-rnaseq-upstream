#!/usr/bin/env bash
#
# count_study.sh -- one count matrix per study, over all of its batches.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# USAGE
#   sbatch count_study.sh <study>
#   sbatch --export=ALL,PAX6_ASSEMBLY=GRCm39 count_study.sh <study>
#
#   The matrix is built from <mapping_root>/<assembly>/ only, and the output
#   directory carries the assembly in its name, so two assemblies produce two
#   matrices and never one mixed one.
#
# THE TWO RULES THIS SCRIPT ENFORCES
#   1. Every library in a study is counted in a SINGLE featureCounts
#      invocation, one annotation, one flag set, one Subread version.
#      Counting previously lived in the per-batch alignment script, so the
#      2025 and 2026 batches were counted months apart in separate jobs, with
#      nothing recording what ran -- a processing difference that could sit
#      exactly on top of the batch covariate used in every model. Adding a batch means re-counting the
#      whole study, which took about 40 minutes for 48 libraries.
#
#   2. Studies are counted SEPARATELY. Lauderdale and GSE183742 carry
#      different Pax6 alleles (Sey-Neu vs tm1Pgr). Quantifying them
#      identically makes them comparable; counting them into one matrix would
#      imply a joint model the design does not support.
#
#SBATCH --job-name=count_study
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=count_%j.out
#SBATCH --error=count_%j.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$PIPELINE_DIR/pipeline_common.sh"

STUDY="${1:-}"
[[ -n "$STUDY" ]] || die "Usage: sbatch count_study.sh <study>"
require_manifest; require_study "$STUDY"

OUT_DIR="${PAX6_COUNTS_OUT:-$PROJECT_ROOT/counts/${STUDY}_${ASSEMBLY}_$(date +%Y%m%d)}"
COUNTS="$OUT_DIR/gene_counts_featureCounts.txt"
BAM_LIST="$OUT_DIR/bam_files.txt"
SAMPLE_MAP="$OUT_DIR/sample_map.tsv"

# --- strand and layout must be measured, agreed and identical -------------
STRAND=""; LAYOUT=""
while read -r B; do
  S=$(manifest_field "$STUDY" "$B" 4)
  L=$(manifest_field "$STUDY" "$B" 3)
  case "$S" in
    0|1|2) ;;
    auto)  die "$STUDY/$B has strand=auto. Run: sbatch strand_check.sh $STUDY $B, then record it in $MANIFEST" ;;
    *)     die "$STUDY/$B has an invalid strand '$S'. Use 0, 1, 2 or auto." ;;
  esac
  case "$L" in
    PE|SE) ;;
    auto)  die "$STUDY/$B has layout=auto. strand_check.sh reports it; record it in $MANIFEST" ;;
    *)     die "$STUDY/$B has an invalid layout '$L'. Use PE, SE or auto." ;;
  esac
  [[ -z "$STRAND" ]] && STRAND="$S" || [[ "$S" == "$STRAND" ]] || die \
"Batches of '$STUDY' disagree on strandedness ($STRAND vs $S).
        One featureCounts run applies one -s to every library. Do not pick one
        and proceed: count the groups separately and handle the difference
        explicitly in the analysis."
  [[ -z "$LAYOUT" ]] && LAYOUT="$L" || [[ "$L" == "$LAYOUT" ]] || die \
"Batches of '$STUDY' mix PE and SE ($LAYOUT vs $L).
        Fragment counts and read counts are not on the same scale and must not
        share a matrix. Split the study, or count the SE batches separately."
done < <(manifest_batches "$STUDY")
log "$STUDY: assembly $ASSEMBLY ($INDEX_TYPE index), layout $LAYOUT, strand -s $STRAND"
log "  annotation: $ANNOTATION"

load_pinned Subread SAMtools
report_versions
[[ -f "$ANNOTATION" ]] || die "Annotation not found: $ANNOTATION"
mkdir -p "$OUT_DIR"

# --- collect and verify every BAM ----------------------------------------
: > "$BAM_LIST"; : > "$SAMPLE_MAP.tmp"; n_total=0
while read -r B; do
  MAPPING_DIR=$(mapping_dir_for "$STUDY" "$B")
  [[ -d "$MAPPING_DIR" ]] || die "$STUDY/$B: no BAMs for assembly $ASSEMBLY at $MAPPING_DIR
        Run: PAX6_ASSEMBLY=$ASSEMBLY sbatch align_batch.sh $STUDY $B"
  n_batch=0
  while read -r BAM; do
    bam_ok "$BAM" || die "$STUDY/$B: BAM incomplete, unindexed or truncated: $BAM
        Re-run: sbatch align_batch.sh $STUDY $B"
    printf '%s\n' "$BAM" >> "$BAM_LIST"
    printf '%s\t%s\t%s\t%s\t%s\n' "$BAM" "$(sample_id_from_bam "$BAM")" "$STUDY" "$B" "$ASSEMBLY" >> "$SAMPLE_MAP.tmp"
    n_batch=$((n_batch+1))
  done < <(find -L "$MAPPING_DIR" -maxdepth 1 -type f -name "*_hisat2.sorted.bam" | sort)
  (( n_batch > 0 )) || die "$STUDY/$B: no BAMs in $MAPPING_DIR"
  log "batch $B: $n_batch libraries"
  n_total=$((n_total+n_batch))
done < <(manifest_batches "$STUDY")

{ printf 'bam_path\tsample_id\tstudy\tbatch\tassembly\n'; cat "$SAMPLE_MAP.tmp"; } > "$SAMPLE_MAP"
rm -f "$SAMPLE_MAP.tmp"
log "$n_total libraries in study $STUDY"

# A sample ID repeated within a study would give two columns that the R loader
# maps to one metadata row. Catch it here, where the fix is a rename.
DUPES=$(awk -F'\t' 'NR>1 { print $2 }' "$SAMPLE_MAP" | sort | uniq -d)
if [[ -n "$DUPES" ]]; then
  echo "[FATAL] Sample IDs occur in more than one batch of '$STUDY':" >&2
  while read -r D; do
    printf '        %s: ' "$D" >&2
    awk -F'\t' -v d="$D" 'NR>1 && $2 == d { printf "%s ", $4 }' "$SAMPLE_MAP" >&2; echo >&2
  done <<< "$DUPES"
  die "Sample IDs must be unique within a study."
fi

[[ -n "${PAX6_EXPECT_N:-}" && "$n_total" != "$PAX6_EXPECT_N" ]] &&
  die "Expected $PAX6_EXPECT_N libraries, found $n_total."
[[ -e "$COUNTS" ]] && die "$COUNTS exists. Move or rename it; this script will not overwrite a count matrix."

# --- count ---------------------------------------------------------------
# --countReadPairs is REQUIRED under Subread >= 2.0.2 to count fragments; -p
# alone only declares paired-end input, and -B/-C then have nothing to act on.
# That is the defect that produced the superseded matrices.
if [[ "$LAYOUT" == PE ]]; then PE_ARGS=(-p --countReadPairs -B -C); else PE_ARGS=(); fi
log "counting $n_total libraries in one invocation"
featureCounts -T "$THREADS" -a "$ANNOTATION" -F GTF \
  "${PE_ARGS[@]}" -t exon -g gene_id \
  -s "$STRAND" -O -M --primary -Q 10 \
  -o "$COUNTS" $(cat "$BAM_LIST")
[[ -s "$COUNTS" ]] || die "featureCounts produced no output"

# --- verify before anyone builds on it -----------------------------------
echo
echo "==================== ASSIGNMENT SUMMARY ===================="
awk -F'\t' '
  NR == 1 { for (i = 2; i <= NF; i++) { n = split($i, p, "/"); name[i] = p[n] } ; next }
  { for (i = 2; i <= NF; i++) { tot[i] += $i; v[$1, i] = $i } }
  END {
    printf "%-46s %10s %12s %10s\n", "Sample", "Assigned%", "NoFeatures%", "A/(A+NF)"
    for (i = 2; i <= NF; i++) {
      a = v["Assigned", i]; nf = v["Unassigned_NoFeatures", i]; r = 100*a/(a+nf)
      printf "%-46s %9.1f%% %11.1f%% %9.1f%%%s\n", name[i], 100*a/tot[i], 100*nf/tot[i], r, (r < 70 ? "   <-- LOW" : "")
      if (r < 70) low++
    }
    if (low) {
      printf "\n%d librar(y/ies) below 70%% on A/(A+NF).\n", low
      print "Near 50% across the board means the strand setting is wrong."
      print "A single low library is more likely degraded or low-complexity RNA."
      exit 3
    }
  }' "${COUNTS}.summary" || RATIO_FAIL=$?
echo "==========================================================="

log "counts     : $COUNTS"
log "sample map : $SAMPLE_MAP"
log "settings   : Subread pinned, $ASSEMBLY, $LAYOUT, -s $STRAND, $n_total libraries"
[[ -n "${RATIO_FAIL:-}" ]] && { warn "Assignment ratios flagged. Resolve before using this matrix."; exit 3; }

cat <<EOF

NEXT
  1. Point the analysis config's FILE_COUNTS at:
       $COUNTS
  2. Archive any superseded matrix OUTSIDE the directory the config searches.
  3. Ensure every sample_id in $SAMPLE_MAP has a metadata row,
     including its batch label.
  4. Re-run the analysis into a NEW output subdirectory, set INSIDE R:
       Sys.setenv(PAX6_OUT_SUBDIR = "ThreeState_$(date +%Y%m%d)"); source("run_all.R")
  5. Diff against the previous run before touching the manuscript.
EOF
