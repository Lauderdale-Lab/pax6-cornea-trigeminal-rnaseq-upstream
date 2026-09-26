#!/usr/bin/env bash
#
# slurm/qc_trim.sbatch -- one library: FastQC (raw), adapter detection,
# Trimmomatic, FastQC (trimmed, both mates). Submitted as an array by
# `bin/pax6 submit`; task N handles line N of the dataset's samples.tsv.
#
#SBATCH --job-name=pax6_trim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -euo pipefail
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
load_run "$PAX6_RUN"
DATASET="$PAX6_DATASET"
IFS=$'\t' read -r SAMPLE R1 R2 < <(sample_line "$(samples_tsv "$DATASET")" "${SLURM_ARRAY_TASK_ID:?not an array task}")

TRIM=$(trim_dir "$DATASET"); QC=$(qc_dir "$DATASET")
mkdir -p "$TRIM" "$QC/raw" "$QC/trimmed" "$QC/adapters" "$QC/trimmomatic"
DONE="$TRIM/$SAMPLE.done"
log "$RUN / $DATASET / $SAMPLE"

if [[ -f "$DONE" ]]; then
  if [[ "$R2" == - ]] && gz_ok "$TRIM/${SAMPLE}.fq.gz"; then log "already trimmed and verified"; exit 0; fi
  if [[ "$R2" != - ]] && gz_ok "$TRIM/${SAMPLE}_1.fq.gz" && gz_ok "$TRIM/${SAMPLE}_2.fq.gz"; then log "already trimmed and verified"; exit 0; fi
  warn "marker present but outputs incomplete; redoing"; rm -f "$DONE"
fi

INPUTS=("$R1"); [[ "$R2" != - ]] && INPUTS+=("$R2")

# --- FastQC on raw reads -------------------------------------------------
load_modules FastQC
fastqc --version
fastqc -q -t "$THREADS" -o "$QC/raw" "${INPUTS[@]}"

# --- Adapter detection ----------------------------------------------------
# Adapter-like overrepresented sequences from FastQC. When none are found the
# TruSeq3 set is used, and that outcome is recorded rather than left silent.
ADAPT_FA="$QC/adapters/${SAMPLE}_detected.fa"; : > "$ADAPT_FA"
for f in "${INPUTS[@]}"; do
  z="$QC/raw/$(basename "${f%.gz}")"; z="${z%.fastq}"; z="${z%.fq}_fastqc.zip"
  [[ -s "$z" ]] || die "FastQC output missing: $z"
  unzip -p "$z" '*/fastqc_data.txt' | awk -F'\t' '
    /^>>Overrepresented sequences/ { m = 1; next }
    m && /^>>END_MODULE/           { exit }
    m && !/^#/ && length($1) >= 12 && $1 !~ /^(A+|C+|G+|T+)$/ &&
      tolower($4) ~ /adapter|illumina|nextera|truseq|nebnext|smallrna/ { print $1 }'
done | sort -u | awk '{ printf ">detected_%d\n%s\n", NR, $0 }' > "$ADAPT_FA"

load_modules Trimmomatic
JAR=$(find "$EBROOTTRIMMOMATIC" -maxdepth 1 -name 'trimmomatic*.jar' | head -n1)
[[ -s "$JAR" ]] || die "Trimmomatic jar not found under $EBROOTTRIMMOMATIC"
if [[ -s "$ADAPT_FA" ]]; then
  ADAPTERS="$ADAPT_FA"; SOURCE=detected
else
  ADAPTERS="$EBROOTTRIMMOMATIC/adapters/TruSeq3-$([[ "$R2" == - ]] && echo SE || echo PE).fa"; SOURCE=TruSeq3_default
fi
printf '%s\t%s\t%s\n' "$SAMPLE" "$SOURCE" "$ADAPTERS" > "$QC/adapters/${SAMPLE}_adapter_source.tsv"
log "adapters: $SOURCE"

# --- Trimmomatic -----------------------------------------------------------
# Written to temporary names and renamed only on success, so a killed job can
# never leave a plausible-looking trimmed file behind. Unpaired survivors are
# not used downstream and are discarded.
TLOG="$QC/trimmomatic/${SAMPLE}.trimmomatic.log"
# shellcheck disable=SC2086  # TRIM_STEPS is a deliberate word list
if [[ "$R2" != - ]]; then
  O1="$TRIM/${SAMPLE}_1.fq.gz"; O2="$TRIM/${SAMPLE}_2.fq.gz"
  java -jar "$JAR" PE -threads "$THREADS" -phred33 "$R1" "$R2" \
    "$O1.tmp.gz" "$TRIM/${SAMPLE}_1.unpaired.tmp.gz" "$O2.tmp.gz" "$TRIM/${SAMPLE}_2.unpaired.tmp.gz" \
    ILLUMINACLIP:"$ADAPTERS":"$TRIM_ILLUMINACLIP" $TRIM_STEPS 2> "$TLOG"
  rm -f "$TRIM/${SAMPLE}"_[12].unpaired.tmp.gz
  { gz_ok "$O1.tmp.gz" && gz_ok "$O2.tmp.gz"; } || die "Trimmomatic output failed verification"
  mv "$O1.tmp.gz" "$O1"; mv "$O2.tmp.gz" "$O2"; OUTS=("$O1" "$O2")
else
  O="$TRIM/${SAMPLE}.fq.gz"
  java -jar "$JAR" SE -threads "$THREADS" -phred33 "$R1" "$O.tmp.gz" \
    ILLUMINACLIP:"$ADAPTERS":"$TRIM_ILLUMINACLIP" $TRIM_STEPS 2> "$TLOG"
  gz_ok "$O.tmp.gz" || die "Trimmomatic output failed verification"
  mv "$O.tmp.gz" "$O"; OUTS=("$O")
fi
grep -E '^(Input Read Pairs|Input Reads):' "$TLOG" || true

# --- FastQC on trimmed reads (every mate) ---------------------------------
load_modules FastQC
fastqc -q -t "$THREADS" -o "$QC/trimmed" "${OUTS[@]}"

touch "$DONE"
log "done"
