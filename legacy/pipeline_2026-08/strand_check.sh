#!/usr/bin/env bash
#
# strand_check.sh -- measure strandedness AND layout for one batch.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# USAGE
#   sbatch strand_check.sh <study> <batch> [n_samples]     # default 3
#
# WHY THIS IS A REQUIRED STEP AND NOT AN ASSUMPTION
#   Both Lauderdale batches were counted for months at -s 2 on the assumption
#   of a directional dUTP library. They are unstranded. Roughly half of every
#   library's exon-overlapping fragments were discarded. The error was
#   invisible downstream precisely because it applied to every sample equally:
#   it distorted no comparison, it only destroyed power.
#
#   For public data the risk is higher, not lower. A GEO record's library
#   description is not evidence, and neither is the vendor.
#
#SBATCH --job-name=strand_check
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=strand_%j.out
#SBATCH --error=strand_%j.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$PIPELINE_DIR/pipeline_common.sh"

STUDY="${1:-}"; BATCH="${2:-}"; N_TEST="${3:-3}"
[[ -n "$STUDY" && -n "$BATCH" ]] || die "Usage: sbatch strand_check.sh <study> <batch> [n]"

require_manifest; require_study "$STUDY"
MAPPING_DIR=$(mapping_dir_for "$STUDY" "$BATCH")
[[ -n "$MAPPING_DIR" ]] || die "No manifest row for $STUDY/$BATCH"
[[ -d "$MAPPING_DIR" ]] || die "Mapping directory not found: $MAPPING_DIR
        (assembly $ASSEMBLY -- run align_batch.sh for this assembly first)"

load_pinned Subread SAMtools
report_versions

OUT_DIR="$MAPPING_DIR/StrandCheck"; mkdir -p "$OUT_DIR"
mapfile -t BAMS < <(find -L "$MAPPING_DIR" -maxdepth 1 -type f -name "*_hisat2.sorted.bam" | sort | head -n "$N_TEST")
(( ${#BAMS[@]} > 0 )) || die "No BAMs in $MAPPING_DIR -- run align_batch.sh $STUDY $BATCH first"

# Layout, measured from the alignment flags.
LAYOUT=$(detect_layout "${BAMS[0]}") || die "Could not read ${BAMS[0]}"
log "layout: $LAYOUT (measured from the paired flag)"
if [[ "$LAYOUT" == PE ]]; then PE_ARGS=(-p --countReadPairs -B -C); else PE_ARGS=(); fi

for BAM in "${BAMS[@]}"; do
  bam_ok "$BAM" || die "BAM incomplete or unindexed: $BAM"
  S=$(sample_id_from_bam "$BAM")
  for X in 0 1 2; do
    OUT="$OUT_DIR/${S}_s${X}.txt"
    [[ -s "${OUT}.summary" ]] && { log "$S -s $X: cached"; continue; }
    log "$S -s $X"
    featureCounts -T "$THREADS" -a "$ANNOTATION" -F GTF \
      "${PE_ARGS[@]}" -t exon -g gene_id \
      -s "$X" -O -M --primary -Q 10 \
      -o "$OUT" "$BAM" > /dev/null 2>&1
  done
done

echo
echo "===== STRAND CHECK: $STUDY / $BATCH  $ASSEMBLY  ($LAYOUT) ====="
printf '%-28s %4s %11s %13s %11s\n' Sample -s "Assigned%" "NoFeatures%" "A/(A+NF)"
for BAM in "${BAMS[@]}"; do
  S=$(sample_id_from_bam "$BAM")
  for X in 0 1 2; do
    awk -v x="$X" -v smp="$S" '
      NR > 1 { t += $2; v[$1] = $2 }
      END { a = v["Assigned"]; nf = v["Unassigned_NoFeatures"]
            printf "%-28s %4s %10.1f%% %12.1f%% %10.1f%%\n", smp, x, 100*a/t, 100*nf/t, 100*a/(a+nf) }' \
      "$OUT_DIR/${S}_s${X}.txt.summary"
  done
done
echo "================================================================"

CALL=$(for BAM in "${BAMS[@]}"; do
  S=$(sample_id_from_bam "$BAM")
  for X in 0 1 2; do
    awk -v x="$X" 'NR>1 { t += $2; v[$1] = $2 } END { printf "%s %.4f\n", x, v["Assigned"]/t }' \
      "$OUT_DIR/${S}_s${X}.txt.summary"
  done
done | awk '
  { sum[$1] += $2; n[$1]++ }
  END {
    a0 = sum[0]/n[0]; a1 = sum[1]/n[1]; a2 = sum[2]/n[2]
    if      (a2 > 0.6*a0 && a1 < 0.25*a0) print "2"
    else if (a1 > 0.6*a0 && a2 < 0.25*a0) print "1"
    else if (a1 > 0.3*a0 && a2 > 0.3*a0 && a1 < 0.7*a0 && a2 < 0.7*a0) print "0"
    else print "AMBIGUOUS"
  }')

echo
if [[ "$CALL" == "AMBIGUOUS" ]]; then
  echo "RESULT: AMBIGUOUS -- the three settings match no expected pattern."
  echo "        Inspect the table and the .summary files in $OUT_DIR."
  echo "        Do not guess."
else
  echo "RESULT: $STUDY/$BATCH is layout = $LAYOUT, strand = $CALL"
  echo "        Set both columns for that row in:"
  echo "        $MANIFEST"
fi
echo
echo "Full outputs: $OUT_DIR"
