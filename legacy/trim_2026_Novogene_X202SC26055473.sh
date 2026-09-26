#!/bin/bash

# Auto-discover samples (nested), support .fastq.gz + .fq.gz and uncompressed .fastq + .fq,
# autodetect adapters from FastQC Overrepresented sequences, then trim + re-QC + MultiQC.
# Prefers newest module versions at runtime.

#SBATCH --job-name=trimmomatic_fastqc_multiqc_adapterauto_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=40:00:00
#SBATCH --output=fastqc_trimmomatic_fastqc_multiqc_adapterauto_all_%j.out
#SBATCH --error=fastqc_trimmomatic_fastqc_multiqc_adapterauto_all_%j.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=BEGIN,FAIL,END


# ---------------- Auto-pick latest module versions ----------------
latest_ver() {
  local family="$1"
  module -t avail "$family" 2>&1 \
    | awk -F'/' '/^'"$family"'\// {print $2}' \
    | sort -V \
    | tail -n1
}

# Get latest versions
TRIMMOD_VER=$(latest_ver Trimmomatic)
FASTQC_VER=$(latest_ver FastQC)
MULTIQC_VER=$(latest_ver MultiQC)

# Load modules (Lmod will handle Java dependency for Trimmomatic)
if [[ -n "$TRIMMOD_VER" ]]; then
  module load "Trimmomatic/$TRIMMOD_VER"
else
  echo "ERROR: No Trimmomatic module found. Exiting."
  exit 1
fi



[[ -n "$FASTQC_VER" ]] && module load "FastQC/$FASTQC_VER" || module load FastQC/0.11.9-Java-17
[[ -n "$MULTIQC_VER" ]] && module load "MultiQC/$MULTIQC_VER" || module load MultiQC/1.14-foss-2022a


# ---------------- Directories ----------------
BASE_DIR="/work/jdllab/Mouse_Pax6_Cornea_2_2026/Novogene_X202SC26055473/01.RawData"
ALL_SAMPLES_LIST="$BASE_DIR/all_discovered_samples.txt"
NEEDS_FASTQC_LIST="$BASE_DIR/needs_fastqc_processing.txt"
RAW_DIR="$BASE_DIR" # supports nested structure; 25 Jun 2026, edited from RAW_DIR="$BASE_DIR/FASTQ_Files"
TRIMMED_DIR="$BASE_DIR/Trimmed_Data"
FASTQC_RAW="$BASE_DIR/FastQC_Raw"
FASTQC_TRIMMED="$BASE_DIR/FastQC_Trimmed"
DETECTED_ADAPTERS="$BASE_DIR/Detected_Adapters"

DEFAULT_ADAPTERS="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
ALT_ADAPTERS_NEXTERA="$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa"
ALT_ADAPTERS_TRUSEQ2="$EBROOTTRIMMOMATIC/adapters/TruSeq2-PE.fa"

LOG_FILE="$BASE_DIR/fastqc_trimmomatic_fastqc_multiqc_adapterauto_run.log"
SUMMARY_REPORT="$BASE_DIR/fastqc_multiqc_summary_adapterauto.txt"
MULTIQC_REPORT="$BASE_DIR/multiqc_report.html"

SKIPPED_FASTQC_RAW="$FASTQC_RAW/skipped_files.log"
SKIPPED_TRIMMOMATIC="$BASE_DIR/skipped_trimmomatic.log"
SKIPPED_FASTQC_TRIMMED="$FASTQC_TRIMMED/skipped_files.log"

RUNTIME_CSV="$BASE_DIR/sample_runtimes.csv"

mkdir -p "$TRIMMED_DIR" "$FASTQC_RAW" "$FASTQC_TRIMMED" "$DETECTED_ADAPTERS"

# Initialize logs
> "$LOG_FILE"; > "$SUMMARY_REPORT"; > "$ALL_SAMPLES_LIST"; > "$NEEDS_FASTQC_LIST"
> "$SKIPPED_FASTQC_RAW"; > "$SKIPPED_TRIMMOMATIC"; > "$SKIPPED_FASTQC_TRIMMED"

echo "Sample,Trimmomatic_Runtime_s,FastQC_Runtime_s,AdapterSource" > "$RUNTIME_CSV"

# Bash 4+ associative arrays
declare -A fastqc_runtimes
declare -A trimmomatic_runtimes
declare -A sample_to_r1
declare -A sample_to_r2
declare -A adapter_source_label

pipeline_start=$(date +%s)
echo "Job started at $(date)" | tee -a "$LOG_FILE"
echo "Scanning RAW_DIR (nested; .fastq.gz/.fq.gz/.fastq/.fq supported): $RAW_DIR" | tee -a "$LOG_FILE"

# ---------------- Discover samples (nested + all common extensions) ----------------
raw_samples=()
declare -A seen_samples
DUPLICATE_LOG="$BASE_DIR/duplicate_samples.log"; > "$DUPLICATE_LOG"; duplicate_found=false

# Log discovered filenames (for visibility only)
find "$RAW_DIR" -type f \
  \( -name "*_1.fastq.gz" -o -name "*_2.fastq.gz" -o -name "*_1.fq.gz" -o -name "*_2.fq.gz" \
   -o -name "*_1.fastq"   -o -name "*_2.fastq"   -o -name "*_1.fq"    -o -name "*_2.fq" \) \
  -printf "%f\n" \
  | sort \
  | tee -a "$LOG_FILE" >/dev/null

# Build sample list from full paths of *_1.* (gz or plain)
while IFS= read -r r1_file; do
  [[ ! -f "$r1_file" ]] && continue
  base_r1=$(basename "$r1_file")

  # Derive sample_name by stripping any of the supported suffixes
  if [[ "$base_r1" == *_1.fastq.gz ]]; then       sample_name="${base_r1%_1.fastq.gz}"
  elif [[ "$base_r1" == *_1.fq.gz ]]; then        sample_name="${base_r1%_1.fq.gz}"
  elif [[ "$base_r1" == *_1.fastq ]]; then        sample_name="${base_r1%_1.fastq}"
  elif [[ "$base_r1" == *_1.fq ]]; then           sample_name="${base_r1%_1.fq}"
  else continue; fi

  r1_dir=$(dirname "$r1_file")
  # Prefer mate in same folder; check all extensions
  mates=(
    "$r1_dir/${sample_name}_2.fastq.gz"
    "$r1_dir/${sample_name}_2.fq.gz"
    "$r1_dir/${sample_name}_2.fastq"
    "$r1_dir/${sample_name}_2.fq"
  )

  if [[ -n "${seen_samples[$sample_name]}" ]]; then
    echo "Duplicate sample name found: $sample_name (existing: ${seen_samples[$sample_name]}, new: $r1_file)" \
      | tee -a "$LOG_FILE" >> "$DUPLICATE_LOG"
    duplicate_found=true
    continue
  fi

  chosen_r2=""
  for cand in "${mates[@]}"; do
    if [[ -f "$cand" ]]; then chosen_r2="$cand"; break; fi
  done

  if [[ -n "$chosen_r2" ]]; then
    seen_samples["$sample_name"]="$r1_file"
    sample_to_r1["$sample_name"]="$r1_file"
    sample_to_r2["$sample_name"]="$chosen_r2"
    raw_samples+=("$sample_name")
  else
    echo "⚠️ Missing mate pair for $sample_name (no _2 file in $r1_dir). Skipping." | tee -a "$LOG_FILE"
  fi
done < <(find "$RAW_DIR" -type f \
          \( -name "*_1.fastq.gz" -o -name "*_1.fq.gz" -o -name "*_1.fastq" -o -name "*_1.fq" \) \
          -print)

if [[ ${#raw_samples[@]} -eq 0 ]]; then
  echo "❌ No paired-end FASTQ files found in $RAW_DIR. Exiting." | tee -a "$LOG_FILE"
  exit 1
fi

! $duplicate_found && echo "No duplicate sample names found." | tee -a "$LOG_FILE"

echo "Discovered ${#raw_samples[@]} samples:" | tee -a "$LOG_FILE"
printf '%s\n' "${raw_samples[@]}" | tee -a "$LOG_FILE" >/dev/null
printf '%s\n' "${raw_samples[@]}" > "$ALL_SAMPLES_LIST"

# ---------------- FastQC on RAW ----------------
echo "Starting FastQC on raw reads..." | tee -a "$LOG_FILE"
missing_samples=()
processed_raw_fastqc=0

for sample_name in "${raw_samples[@]}"; do
  r1_file="${sample_to_r1[$sample_name]}"
  r2_file="${sample_to_r2[$sample_name]}"

  output_r1_zip="$FASTQC_RAW/${sample_name}_1_fastqc.zip"
  output_r2_zip="$FASTQC_RAW/${sample_name}_2_fastqc.zip"

  if [[ ! -s "$output_r1_zip" || ! -s "$output_r2_zip" ]]; then
    missing_samples+=("$sample_name")
    echo "$r1_file" >> "$NEEDS_FASTQC_LIST"
  fi
done

for sample_name in "${missing_samples[@]}"; do
  r1_file="${sample_to_r1[$sample_name]}"
  r2_file="${sample_to_r2[$sample_name]}"

  if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
    echo "❌ Missing inputs for $sample_name" | tee -a "$LOG_FILE"
    continue
  fi

  echo "Processing FastQC for $sample_name..." | tee -a "$LOG_FILE"
  start_time=$(date +%s)
  fastqc -t 4 -o "$FASTQC_RAW" "$r1_file" "$r2_file" >> "$LOG_FILE" 2>&1
  end_time=$(date +%s); elapsed=$((end_time - start_time))
  fastqc_runtimes["$sample_name"]=$elapsed

  if [[ -s "$FASTQC_RAW/${sample_name}_1_fastqc.zip" && -s "$FASTQC_RAW/${sample_name}_2_fastqc.zip" \
        && -s "$FASTQC_RAW/${sample_name}_1_fastqc.html" && -s "$FASTQC_RAW/${sample_name}_2_fastqc.html" ]]; then
    echo "FastQC for $sample_name completed in ${elapsed}s" | tee -a "$LOG_FILE"
    ((processed_raw_fastqc++))
  else
    echo "FastQC failed or incomplete for $sample_name (runtime: ${elapsed}s)" | tee -a "$LOG_FILE"
  fi
done

total_samples=${#raw_samples[@]}
skipped_raw_fastqc=$(( total_samples - processed_raw_fastqc ))

# ---------------- Adapter autodetection from FastQC ----------------
echo "Autodetecting adapters from FastQC Overrepresented sequences..." | tee -a "$LOG_FILE"

extract_overrep() {
  local zipfile="$1"
  # Pull fastqc_data.txt and emit "sequence<TAB>source" lines for the Overrepresented sequences module
  unzip -p "$zipfile" */fastqc_data.txt 2>/dev/null \
  | awk '
    BEGIN{ inmod=0 }
    /^>>Overrepresented sequences/ { inmod=1; next }
    inmod==1 && /^>>END_MODULE/   { inmod=0; exit }
    inmod==1 && $0!~/^#/ {
      n=split($0,a,"\t");
      if (n>=4) { seq=a[1]; src=a[4]; print seq"\t"src }
    }'
}

for sample_name in "${raw_samples[@]}"; do
  r1_zip="$FASTQC_RAW/${sample_name}_1_fastqc.zip"
  r2_zip="$FASTQC_RAW/${sample_name}_2_fastqc.zip"
  tmp_overrep="$DETECTED_ADAPTERS/${sample_name}_overrep.tsv"; > "$tmp_overrep"

  [[ -s "$r1_zip" ]] && extract_overrep "$r1_zip" >> "$tmp_overrep"
  [[ -s "$r2_zip" ]] && extract_overrep "$r2_zip" >> "$tmp_overrep"

  cleaned_tsv="$DETECTED_ADAPTERS/${sample_name}_overrep_clean.tsv"
  awk 'BEGIN{IGNORECASE=1}
       {
         seq=toupper($1); src=$2;
         # Keep plausible adapter-like sequences (>=12nt) and drop homopolymers
         if (length(seq)>=12 && seq !~ /^(A{12,}|T{12,}|C{12,}|G{12,})$/) {
           if (src ~ /(Adapter|Illumina|Nextera|TruSeq|NEBNext|SmallRNA)/) {
             print seq"\t"src
           }
         }
       }' "$tmp_overrep" \
    | sort -u > "$cleaned_tsv"

  adapters_fa="$DETECTED_ADAPTERS/${sample_name}_adapters.fa"; > "$adapters_fa"; idx=0
  while IFS=$'\t' read -r seq src; do
    [[ -z "$seq" ]] && continue
    idx=$((idx+1))
    hdr_src=$(echo "$src" | tr ' ' '_' | tr -cd '[:alnum:]_')
    echo ">${sample_name}.adapter_${idx}.${hdr_src}" >> "$adapters_fa"
    echo "$seq" >> "$adapters_fa"
  done < "$cleaned_tsv"

  if [[ -s "$adapters_fa" ]]; then
    adapter_source_label["$sample_name"]="Detected"
    echo "[$sample_name] Using detected adapters: $adapters_fa" | tee -a "$LOG_FILE"
  else
    adapter_source_label["$sample_name"]="Default_TruSeq3-PE"
    echo "[$sample_name] No adapters detected; using default: $DEFAULT_ADAPTERS" | tee -a "$LOG_FILE"
  fi
done

# ---------------- Trimmomatic (per-sample adapters) ----------------
echo "Running Trimmomatic with per-sample adapters..." | tee -a "$LOG_FILE"
processed_trim=0; skipped_trim=0

for sample_name in "${raw_samples[@]}"; do
  R1="${sample_to_r1[$sample_name]}"; R2="${sample_to_r2[$sample_name]}"
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "❌ Missing pair for $sample_name" | tee -a "$LOG_FILE"
    continue
  fi

  out_prefix="$TRIMMED_DIR/$sample_name"
  out1="${out_prefix}_1_paired.fq.gz"
  out2="${out_prefix}_2_paired.fq.gz"
  out1u="${out_prefix}_1_unpaired.fq.gz"
  out2u="${out_prefix}_2_unpaired.fq.gz"

  if [[ -s "$out1" && -s "$out2" && -s "$out1u" && -s "$out2u" ]]; then
    echo "Skipping $sample_name — trimmed outputs already exist." | tee -a "$LOG_FILE"
    echo "$R1" >> "$SKIPPED_TRIMMOMATIC"; echo "$R2" >> "$SKIPPED_TRIMMOMATIC"
    ((skipped_trim++))
    continue
  fi

  adapters_fa="$DETECTED_ADAPTERS/${sample_name}_adapters.fa"
  [[ ! -s "$adapters_fa" ]] && adapters_fa="$DEFAULT_ADAPTERS"

  echo "Processing $sample_name with Trimmomatic (adapters: $adapters_fa)..." | tee -a "$LOG_FILE"
  start_time=$(date +%s)
  java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" \
    PE -threads 8 -phred33 \
    "$R1" "$R2" \
    "$out1" "$out1u" "$out2" "$out2u" \
    ILLUMINACLIP:"$adapters_fa":2:30:10:2:True \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
    >> "$LOG_FILE" 2>&1
  end_time=$(date +%s); elapsed=$((end_time - start_time))
  trimmomatic_runtimes["$sample_name"]=$elapsed
  echo "Finished $sample_name in ${elapsed}s" | tee -a "$LOG_FILE"
  ((processed_trim++))
done

# ---------------- FastQC after trimming ----------------
echo "Running FastQC on trimmed reads..." | tee -a "$LOG_FILE"
fastqc_failed=0

# Collect trimmed R1 paired files
mapfile -d '' trimmed_files < <(find "$TRIMMED_DIR" -type f -name "*_1_paired.fq.gz" -print0)
files_to_process=()

for R1 in "${trimmed_files[@]}"; do
  [[ ! -f "$R1" ]] && { echo "❌ Skipping $R1 — file not found or unreadable" | tee -a "$LOG_FILE"; continue; }
  base=$(basename "$R1" .fq.gz)
  report_html="$FASTQC_TRIMMED/${base}_fastqc.html"
  report_zip="$FASTQC_TRIMMED/${base}_fastqc.zip"
  if [[ -s "$report_html" && -s "$report_zip" ]]; then
    echo "Skipping FastQC for $base (report exists)" | tee -a "$LOG_FILE"
    echo "$R1" >> "$SKIPPED_FASTQC_TRIMMED"
  else
    files_to_process+=("$R1")
  fi
done

for R1 in "${files_to_process[@]}"; do
  base=$(basename "$R1" .fq.gz)
  echo "Running FastQC for $base..." | tee -a "$LOG_FILE"
  start_time=$(date +%s)
  fastqc -o "$FASTQC_TRIMMED" "$R1" >> "$LOG_FILE" 2>&1
  end_time=$(date +%s); elapsed=$((end_time - start_time))
  report_html="$FASTQC_TRIMMED/${base}_fastqc.html"
  report_zip="$FASTQC_TRIMMED/${base}_fastqc.zip"
  if [[ -s "$report_html" && -s "$report_zip" ]]; then
    echo "FastQC completed successfully for $base in ${elapsed}s" | tee -a "$LOG_FILE"
  else
    echo "FastQC failed or incomplete for $base (runtime: ${elapsed}s)" | tee -a "$LOG_FILE"
    ((fastqc_failed++))
  fi
done

# ---------------- MultiQC summary ----------------
echo "Running MultiQC..." | tee -a "$LOG_FILE"
multiqc "$BASE_DIR" -o "$BASE_DIR" >> "$LOG_FILE" 2>&1

if [[ -f "$MULTIQC_REPORT" ]]; then
  multiqc_status="Success"
  echo "MultiQC report generated: $MULTIQC_REPORT" | tee -a "$LOG_FILE"
else
  multiqc_status="Failed"
  echo "MultiQC report not found. MultiQC may have failed." | tee -a "$LOG_FILE"
fi

pipeline_end=$(date +%s)
total_runtime=$((pipeline_end - pipeline_start))

# ---------------- Summary report ----------------
{
  echo "Trimming + FastQC + Adapter Autodetect + MultiQC Summary (newest modules; all extensions)"
  echo "-------------------------------------------------------------------------------------------"
  echo "Date: $(date)"
  echo ""
  echo "Directories:"
  echo "  Raw data: $RAW_DIR (nested; .fastq.gz/.fq.gz/.fastq/.fq)"
  echo "  Trimmed data: $TRIMMED_DIR"
  echo "  FastQC (before trimming): $FASTQC_RAW"
  echo "  FastQC (after trimming):  $FASTQC_TRIMMED"
  echo "  Detected adapters:        $DETECTED_ADAPTERS"
  echo "  MultiQC report:           $MULTIQC_REPORT"
  echo "  MultiQC status:           $multiqc_status"
  echo ""
  echo "Additional Reports:"
  echo "  All discovered SRR samples: $ALL_SAMPLES_LIST"
  echo "  Samples needing FastQC:     $NEEDS_FASTQC_LIST"
  echo ""
  echo "Raw FastQC:"
  echo "  Total raw samples found:    $total_samples"
  echo "  Raw samples processed:      $processed_raw_fastqc"
  echo "  Raw samples skipped:        $skipped_raw_fastqc"
  echo "  Skipped raw FastQC files log: $SKIPPED_FASTQC_RAW"
  echo ""
  echo "Trimmomatic:"
  echo "  Trimmed files processed:    $processed_trim"
  echo "  Trimmed files skipped:      $skipped_trim"
  echo "  Skipped Trimmomatic files log: $SKIPPED_TRIMMOMATIC"
  echo ""
  echo "FastQC after trimming:"
  echo "  FastQC failures after trimming: $fastqc_failed"
  echo "  Skipped trimmed FastQC files log: $SKIPPED_FASTQC_TRIMMED"
  echo ""
  echo "Performance:"
  echo "  Total pipeline runtime:     ${total_runtime}s"
  echo "  Per-sample runtimes saved to: $RUNTIME_CSV"
} > "$SUMMARY_REPORT"

# Append adapter source to CSV
for sample in "${!trimmomatic_runtimes[@]}"; do
  echo "$sample,${trimmomatic_runtimes[$sample]},${fastqc_runtimes[$sample]:-NA},${adapter_source_label[$sample]:-Default_TruSeq3-PE}" \
    >> "$RUNTIME_CSV"
done

# Email the summary
(
  echo "Subject: FastQC_raw + AdapterAuto + Trimming + FastQC_Trimmed + MultiQC (Newest Modules; All Extensions)"
  echo "To: jdlauder@uga.edu"
  echo "Content-Type: text/plain"
  echo ""
  cat "$SUMMARY_REPORT"
) | /usr/sbin/sendmail -t

echo "Pipeline completed at $(date)" | tee -a "$LOG_FILE"
