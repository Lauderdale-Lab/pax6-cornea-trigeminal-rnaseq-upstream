#!/usr/bin/env bash
#
# sam_to_bam.sh -- convert archived SAM files to sorted, indexed BAMs.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# USAGE
#   sbatch sam_to_bam.sh <study> <batch> <sam_dir>
#
# READ THIS BEFORE USING THE RESULT
#   Converting a SAM inherits whatever alignment produced it. If those SAMs
#   were made against a different genome build, a different annotation, or
#   different HISAT2 settings than your own libraries, the resulting counts are
#   NOT comparable with yours no matter how carefully they are counted
#   afterwards -- and nothing downstream will warn you.
#
#   This script therefore prints each SAM's @PG line, which records the exact
#   command that produced it, and refuses to continue unless you have looked.
#   Check it yourself first:
#       samtools view -H <file>.sam | grep '^@PG'
#       samtools view -H <file>.sam | grep '^@SQ' | head -3
#
#   If the reference or parameters differ from this project's, re-align from
#   FASTQ instead: fetch the runs from SRA and use qc_trim_batch.sh followed by
#   align_batch.sh. Conversion is for a quick look; re-alignment is for a
#   result you intend to publish.
#
#SBATCH --job-name=sam_to_bam
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=sam2bam_%j.out
#SBATCH --error=sam2bam_%j.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$PIPELINE_DIR/pipeline_common.sh"

STUDY="${1:-}"; BATCH="${2:-}"; SAM_DIR="${3:-}"
[[ -n "$STUDY" && -n "$BATCH" && -n "$SAM_DIR" ]] ||
  die "Usage: sbatch sam_to_bam.sh <study> <batch> <sam_dir>"

require_manifest; require_study "$STUDY"
# Assembly-aware, exactly as align_batch.sh: BAMs must land where
# count_study.sh looks for them, i.e. <mapping_root>/<assembly>/.
[[ -n "$(manifest_field "$STUDY" "$BATCH" 6)" ]] || die "No manifest row for $STUDY/$BATCH"
MAPPING_DIR=$(mapping_dir_for "$STUDY" "$BATCH")
[[ -d "$SAM_DIR" ]]     || die "SAM directory not found: $SAM_DIR"

load_pinned SAMtools
report_versions
mkdir -p "$MAPPING_DIR"

mapfile -t SAMS < <(find -L "$SAM_DIR" -maxdepth 1 -type f \( -name "*.sam" -o -name "*.sam.gz" \) | sort)
(( ${#SAMS[@]} > 0 )) || die "No SAM files in $SAM_DIR"
log "${#SAMS[@]} SAM file(s)"

# --- provenance of the alignment being inherited -------------------------
echo
echo "=============== ALIGNMENT PROVENANCE OF THESE SAMs ==============="
REF_SIG=""
for SAM in "${SAMS[@]}"; do
  echo "--- $(basename "$SAM")"
  samtools view -H "$SAM" 2>/dev/null | grep '^@PG' | head -n2 | sed 's/^/    /' || echo "    (no @PG line)"
  SIG=$(samtools view -H "$SAM" 2>/dev/null | awk '$0 ~ /^@SQ/ { n++; if (n <= 3) printf "%s;", $2 } END { printf "n=%d", n }')
  echo "    reference: $SIG"
  [[ -z "$REF_SIG" ]] && REF_SIG="$SIG"
  [[ "$SIG" == "$REF_SIG" ]] || warn "$(basename "$SAM") has a DIFFERENT reference signature from the first file"
done
echo "=================================================================="
echo
echo "Compare the @PG command and the reference above with this project's"
echo "selected assembly ($ASSEMBLY):"
echo "  index      : $HISAT2_INDEX"
echo "  annotation : $ANNOTATION"
echo "  BAMs go to : $MAPPING_DIR"
echo

if [[ "${PAX6_SAM_PROVENANCE_OK:-}" != "yes" ]]; then
  cat >&2 <<EOF
[STOP] Provenance not confirmed.

  The header above is printed so a person decides whether these alignments
  belong in this project. If the reference and parameters match, re-submit with:

      sbatch --export=ALL,PAX6_SAM_PROVENANCE_OK=yes sam_to_bam.sh $STUDY $BATCH $SAM_DIR

  If they do NOT match, do not convert. Re-align from FASTQ instead.
EOF
  exit 2
fi

n_done=0; n_skip=0
for SAM in "${SAMS[@]}"; do
  S=$(sample_id_from_bam "$SAM")
  BAM="$MAPPING_DIR/${S}_hisat2.sorted.bam"
  if bam_ok "$BAM"; then log "$S: complete BAM present, skipping"; n_skip=$((n_skip+1)); continue; fi
  [[ -e "$BAM" ]] && { warn "$S: incomplete BAM, redoing"; rm -f "$BAM" "$BAM.bai"; }
  log "$S: sorting to BAM"
  samtools sort -@ "$THREADS" -m 2G -o "$BAM" "$SAM"
  samtools index -@ 4 "$BAM"
  bam_ok "$BAM" || die "$S: BAM failed verification"
  n_done=$((n_done+1))
done

log "$n_done converted, $n_skip already present -> $MAPPING_DIR"
log "These BAMs are recorded as assembly $ASSEMBLY. If the SAM headers showed a"
log "different reference, that label is WRONG -- delete them and re-align."
log "The SAMs are not deleted. Remove them yourself once the BAMs are verified."
log "NEXT: sbatch strand_check.sh $STUDY $BATCH   (it also reports PE vs SE)"
