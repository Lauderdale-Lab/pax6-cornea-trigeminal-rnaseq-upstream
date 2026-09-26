#!/usr/bin/env bash
#
# qc_trim_batch.sh -- FastQC, adapter trimming, re-QC and MultiQC for one batch.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# USAGE
#   sbatch qc_trim_batch.sh <study> <batch> <raw_dir>
#
# Writes trimmed FASTQs to the batch's trimmed_dir from studies.tsv, and QC
# reports beside it.
#
# CHANGES FROM THE 2026 TRIMMING SCRIPT
#   * Module versions pinned, not "latest available".
#   * Post-trim FastQC covers BOTH mates. The old script re-QC'd R1 only, so
#     the trimmed-read section of every MultiQC report described half the data.
#   * Single-end input supported; GEO datasets are frequently SE.
#   * Adapter autodetection retained, but its outcome is reported plainly. In
#     the 2026 batch it detected nothing in all 21 libraries and every one fell
#     back to TruSeq3 -- which is a finding worth stating, not a silent default.
#
#SBATCH --job-name=qc_trim
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=qc_trim_%A_%a.out
#SBATCH --error=qc_trim_%A_%a.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$PIPELINE_DIR/pipeline_common.sh"

STUDY="${1:-}"; BATCH="${2:-}"; RAW_DIR="${3:-}"
[[ -n "$STUDY" && -n "$BATCH" && -n "$RAW_DIR" ]] ||
  die "Usage: sbatch qc_trim_batch.sh <study> <batch> <raw_dir>"

require_manifest; require_study "$STUDY"
TRIMMED_DIR=$(manifest_field "$STUDY" "$BATCH" 5)
LAYOUT=$(manifest_field "$STUDY" "$BATCH" 3)
[[ -n "$TRIMMED_DIR" ]] || die "No manifest row for $STUDY/$BATCH"
[[ -d "$RAW_DIR" ]]     || die "Raw directory not found: $RAW_DIR"

load_pinned FastQC Trimmomatic MultiQC
report_versions

QC_RAW="$(dirname "$TRIMMED_DIR")/qc/raw"
QC_TRIM="$(dirname "$TRIMMED_DIR")/qc/trimmed"
ADAPT_DIR="$(dirname "$TRIMMED_DIR")/qc/adapters"
mkdir -p "$TRIMMED_DIR" "$QC_RAW" "$QC_TRIM" "$ADAPT_DIR"

ADAPTERS_PE="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
ADAPTERS_SE="$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE.fa"
TRIM_JAR=$(find "$EBROOTTRIMMOMATIC" -maxdepth 1 -name "trimmomatic*.jar" | head -n1)
[[ -s "$TRIM_JAR" ]] || die "Trimmomatic jar not found under $EBROOTTRIMMOMATIC"

# Layout: measured from the files present, not assumed from the manifest.
# Searched to depth 4 because vendors nest one directory per sample
# (Novogene: 01.RawData/<sample>/<sample>_1.fq.gz).
RAW_DEPTH="${PAX6_RAW_DEPTH:-4}"
mapfile -t R1S < <(find_r1 "$RAW_DIR" "$RAW_DEPTH")
mapfile -t SES < <(find_se "$RAW_DIR" "$RAW_DEPTH")
if (( ${#R1S[@]} > 0 )); then
  FOUND_LAYOUT=PE
elif (( ${#SES[@]} > 0 )); then
  FOUND_LAYOUT=SE
else
  die "No FASTQ files found in $RAW_DIR"
fi
log "$STUDY/$BATCH: $FOUND_LAYOUT input"

# Searching several levels deep means the same sample name could turn up
# twice (a re-delivery, a copy in a subfolder). Two inputs writing one trimmed
# file would be silent corruption, so refuse.
if [[ "$FOUND_LAYOUT" == PE ]]; then _IN=("${R1S[@]}"); else _IN=("${SES[@]}"); fi
DUP_SAMPLES=$(for f in "${_IN[@]}"; do sample_from_r1 "$f"; done | sort | uniq -d)
[[ -z "$DUP_SAMPLES" ]] || die "Sample name(s) found more than once under $RAW_DIR:
$(printf '        %s\n' $DUP_SAMPLES)
        Remove or move the duplicate files before trimming."
log "${#_IN[@]} sample(s)"
if [[ "$LAYOUT" != "auto" && "$LAYOUT" != "$FOUND_LAYOUT" ]]; then
  die "Manifest says layout=$LAYOUT but $RAW_DIR contains $FOUND_LAYOUT data.
        Fix the manifest rather than the data."
fi

# --- FastQC on raw -------------------------------------------------------
log "FastQC on raw reads"
if [[ "$FOUND_LAYOUT" == PE ]]; then
  for R1 in "${R1S[@]}"; do
    R2=$(mate_of "$R1") || die "No mate for $R1"
    fastqc -t "$THREADS" -o "$QC_RAW" "$R1" "$R2"
  done
else
  fastqc -t "$THREADS" -o "$QC_RAW" "${SES[@]}"
fi

# --- Adapter autodetection ----------------------------------------------
# Pull adapter-like overrepresented sequences from the FastQC reports. Report
# the outcome per library: "detected nothing, used the default" is information.
detect_adapters() {
  local sample="$1"; shift
  local out="$ADAPT_DIR/${sample}_adapters.fa"; : > "$out"
  local zip i=0
  for zip in "$@"; do
    [[ -s "$zip" ]] || continue
    unzip -p "$zip" '*/fastqc_data.txt' 2>/dev/null | awk '
      /^>>Overrepresented sequences/ { inmod = 1; next }
      inmod && /^>>END_MODULE/       { exit }
      inmod && $0 !~ /^#/ {
        n = split($0, a, "\t")
        if (n >= 4 && length(a[1]) >= 12 && a[1] !~ /^(A+|T+|C+|G+)$/ &&
            a[4] ~ /(Adapter|Illumina|Nextera|TruSeq|NEBNext|SmallRNA)/)
          print a[1]
      }'
  done | sort -u | while read -r seq; do
    i=$((i + 1)); printf '>%s.detected_%d\n%s\n' "$sample" "$i" "$seq" >> "$out"
  done
  [[ -s "$out" ]] && printf '%s\n' "$out" || printf '\n'
}

# --- Trim ----------------------------------------------------------------
n_detected=0; n_default=0
if [[ "$FOUND_LAYOUT" == PE ]]; then
  for R1 in "${R1S[@]}"; do
    S=$(sample_from_r1 "$R1"); R2=$(mate_of "$R1")
    O1="$TRIMMED_DIR/${S}_1_paired.fq.gz"; O1U="$TRIMMED_DIR/${S}_1_unpaired.fq.gz"
    O2="$TRIMMED_DIR/${S}_2_paired.fq.gz"; O2U="$TRIMMED_DIR/${S}_2_unpaired.fq.gz"
    [[ -s "$O1" && -s "$O2" ]] && { log "$S: trimmed already"; continue; }

    AD=$(detect_adapters "$S" "$QC_RAW/${S}_1_fastqc.zip" "$QC_RAW/${S}_2_fastqc.zip")
    if [[ -n "$AD" ]]; then n_detected=$((n_detected+1)); else AD="$ADAPTERS_PE"; n_default=$((n_default+1)); fi

    log "$S: trimming (adapters: $(basename "$AD"))"
    java -jar "$TRIM_JAR" PE -threads "$THREADS" -phred33 \
      "$R1" "$R2" "$O1" "$O1U" "$O2" "$O2U" \
      ILLUMINACLIP:"$AD":2:30:10:2:True \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
  done
else
  for FQ in "${SES[@]}"; do
    S=$(sample_from_r1 "$FQ")
    OUT="$TRIMMED_DIR/${S}_trimmed.fq.gz"
    [[ -s "$OUT" ]] && { log "$S: trimmed already"; continue; }

    AD=$(detect_adapters "$S" "$QC_RAW/${S}_fastqc.zip")
    if [[ -n "$AD" ]]; then n_detected=$((n_detected+1)); else AD="$ADAPTERS_SE"; n_default=$((n_default+1)); fi

    log "$S: trimming (adapters: $(basename "$AD"))"
    java -jar "$TRIM_JAR" SE -threads "$THREADS" -phred33 "$FQ" "$OUT" \
      ILLUMINACLIP:"$AD":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
  done
fi
log "adapters: $n_detected library(ies) detected, $n_default used the default"

# --- FastQC on trimmed, BOTH mates --------------------------------------
log "FastQC on trimmed reads (both mates)"
mapfile -t TRIMMED < <(find "$TRIMMED_DIR" -maxdepth 1 -type f \
  \( -name "*_paired.fq.gz" -o -name "*_trimmed.fq.gz" \) ! -name "*_unpaired.fq.gz" | sort)
(( ${#TRIMMED[@]} > 0 )) || die "No trimmed FASTQs produced"
fastqc -t "$THREADS" -o "$QC_TRIM" "${TRIMMED[@]}"

# --- MultiQC -------------------------------------------------------------
MQC_DIR="$(dirname "$TRIMMED_DIR")/qc"
log "MultiQC"
multiqc -f -o "$MQC_DIR" -n "multiqc_${STUDY}_${BATCH}.html" "$QC_RAW" "$QC_TRIM"

log "done. QC: $MQC_DIR   trimmed: $TRIMMED_DIR"
log "NEXT: sbatch align_batch.sh $STUDY $BATCH"
