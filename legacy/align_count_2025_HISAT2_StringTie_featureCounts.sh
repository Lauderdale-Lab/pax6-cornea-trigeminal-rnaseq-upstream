#!/usr/bin/env bash



# SLURM job scheduler directives for resource allocation
# Change job name to match project as needed

#SBATCH --job-name=HISAT2_Mouse_RNAseq_Pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=HISAT2_Mouse_RNAseq_%j.out
#SBATCH --error=HISAT2_Mouse_RNAseq_%j.err
#SBATCH --mail-user=jdlauder@uga.edu
#SBATCH --mail-type=ALL

echo -e "\n** Script started on $(date) **\n"

cd "$SLURM_SUBMIT_DIR"

##########################################
# === Error Handling and Logging Setup ===
###########################################
# Enables strict error handling for safer execution:
# -e : Exit immediately if a command exits with a non-zero status.
# -u : Treat unset variables as an error and exit immediately.
# -o pipefail : Return the exit status of the last command in the pipeline that failed.

# Enable strict error handling for safer script execution
set -euo pipefail
set -x

# Redirect all standard error output to a log file for debugging
#Using same logfile as for structured trap messages for simplicity
exec 2>> HISAT2_pipeline_error.log

# Trap and log any command failure with line number and exit code
trap 'echo "[ERROR] Script failed at line $LINENO with exit code $?" >> HISAT2_pipeline_error.log' ERR


####################################
# === USER CONFIGURATION SECTION ===
####################################

# Set the following variables to match your organism and dataset
BASE_DIR_PATH="/scratch/jdlauder/files/Mouse_Cornea_and_Trigeminal_May_2025"  # <-- Update this path to where RNAseq files are and where output is to be written
GENOME_DIR_PATH="/work/jdllab/HISAT2_Mouse_Index/GRCm38"  # <-- Update this path
GENOME_ASSEMBLY="GRCm38.primary_assembly.genome"  # e.g., GRCh38.primary_assembly.genome (without .fa/.fna)
ANNOTATION_FILE="gencode.vM25.annotation.gtf"  # e.g., gencode.v38.annotation.gtf
VCF_FILE="mus_musculus.vcf"  # Optional: SNP VCF file for SNP-aware alignment


###############################################
# === Derived Paths from User Configuration ===
###############################################

GENOME_ASSEMBLY_PATH="$GENOME_DIR_PATH/$GENOME_ASSEMBLY"
ANNOTATION="$GENOME_DIR_PATH/$ANNOTATION_FILE"
VCF_FILE_PATH="$GENOME_DIR_PATH/$VCF_FILE"

RNA_SEQ_TRIMMED_FASTQ_PATH="$BASE_DIR_PATH/Trimmed_Data"
RNA_SEQ_MAPPING_OUTPUT_DIR="$BASE_DIR_PATH/RNAseq_Mapping"

##########################################################
# Dynamically detect HISAT2 index path based on .ht2 files
##########################################################

echo "[INFO] Scanning for valid HISAT2 index sets under $GENOME_DIR_PATH..." | tee -a HISAT2_pipeline_error.log

if [[ "${BASH_VERSINFO:-0}" -lt 4 ]]; then
  echo "[ERROR] This script requires Bash version 4 or higher." >&2
  exit 1
fi

declare -A index_counts
echo "[DEBUG] index_counts associative array declared" | tee -a HISAT2_pipeline_error.log


# Find all .ht2 files and group by full prefix
while IFS= read -r -d '' file; do
  prefix="${file%.*.ht2}"  # Remove .1.ht2, .2.ht2, etc.
  if [[ -n "$prefix" ]]; then
#    ((index_counts["$prefix"]++))
	index_counts["$prefix"]=$(( ${index_counts["$prefix"]:-0} + 1 ))
  else
    echo "[WARNING] Skipping file with empty prefix: $file" | tee -a HISAT2_pipeline_error.log
  fi
done < <(find "$GENOME_DIR_PATH" -type f -name "*.ht2" -print0)


VALID_INDEX_PATHS=()
PREFERRED_INDEX_PATH=""
FOUND_SNP_AWARE=0

for prefix in "${!index_counts[@]}"; do
  echo "[DEBUG] $prefix has ${index_counts[$prefix]} .ht2 files" | tee -a HISAT2_pipeline_error.log
    count=${index_counts[$prefix]}
    if [[ $count -eq 8 ]]; then
        echo "[INFO] Found valid HISAT2 index: $prefix" | tee -a HISAT2_pipeline_error.log
        VALID_INDEX_PATHS+=("$prefix")
        if [[ "$prefix" == *"genome_snp_tran" ]]; then
            echo "[INFO] This index is SNP-aware (genome_snp_tran)." | tee -a HISAT2_pipeline_error.log
            PREFERRED_INDEX_PATH="$prefix"
            FOUND_SNP_AWARE=1
        elif [[ "$prefix" == *"genome_tran" ]]; then
            echo "[INFO] This index is transcript-aware only (genome_tran)." | tee -a HISAT2_pipeline_error.log
        else
            echo "[WARNING] Index not in expected folder (genome_snp_tran or genome_tran)." | tee -a HISAT2_pipeline_error.log
        fi
    else
        echo "[WARNING] Incomplete index set: $prefix ($count files)" | tee -a HISAT2_pipeline_error.log
    fi
done

if [[ ${#VALID_INDEX_PATHS[@]} -eq 0 ]]; then
    echo "[ERROR] No valid HISAT2 index sets (.ht2 files) found under $GENOME_DIR_PATH." | tee -a HISAT2_pipeline_error.log
    exit 1
fi

if [[ $FOUND_SNP_AWARE -eq 0 ]]; then
    PREFERRED_INDEX_PATH="${VALID_INDEX_PATHS[0]}"
fi

if [[ ${#VALID_INDEX_PATHS[@]} -gt 1 ]]; then
    echo "[WARNING] Multiple index sets found. Defaulting to SNP-aware index if available." | tee -a HISAT2_pipeline_error.log
fi

HISAT2_GENOME_INDEX_PATH="$PREFERRED_INDEX_PATH"
HISAT2_GENOME_INDEX_BUILD_LOG="$HISAT2_GENOME_INDEX_PATH/pipeline_index_building.log"
mkdir -p "$(dirname "$HISAT2_GENOME_INDEX_BUILD_LOG")"

echo "[INFO] Selected HISAT2 index path: $HISAT2_GENOME_INDEX_PATH" | tee -a HISAT2_pipeline_error.log

# === Output Directory Setup ===
GTF_RNA_SEQ_MAPPING_OUTPUT_DIR="$BASE_DIR_PATH/StringTie_Assemblies"
COUNTS_OUT="$BASE_DIR_PATH/RNAseq_Quantification_Matrices"

# Create output directories for HISAT2, StringTie, and featureCounts
mkdir -p "$RNA_SEQ_MAPPING_OUTPUT_DIR" "$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR" "$COUNTS_OUT"
touch "$RNA_SEQ_MAPPING_OUTPUT_DIR/test.tmp" || { echo "[ERROR] Cannot write to $RNA_SEQ_MAPPING_OUTPUT_DIR"; exit 1; } | tee -a HISAT2_pipeline_error.log
rm "$RNA_SEQ_MAPPING_OUTPUT_DIR/test.tmp"



###############################################
#Scripting notes:
#when you use \ to have line continuation of a command, there can be no extra tab,
# blank space, or comment after the \

#Defining paths to files
#NOTE: When you define a variable, there can be no blank spaces around the equal sign
#NOTE: as per HISAT2 documentation, do not include trailing ".X.ht2" in filename prefix.
#If included, program will return (ERR): "/(path)/Asag_genome_tran.1.ht2" does not exist
######################################


# Load required software modules for the pipeline
ml purge
ml HISAT2/2.2.1-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0
ml StringTie/2.2.1-GCC-11.3.0
ml Subread/2.0.6-GCC-11.3.0  # for featureCounts

# === Tool Availability Check and Strand Specificity Setting ===
# Ensure all required tools are available in the environment before proceeding.
# Each command checks if the tool is in the system's PATH using `command -v`.
# If any tool is missing, the script exits with an error message.

# Check if all required tools are available in the environment
command -v hisat2 >/dev/null 2>&1 || { echo "[ERROR] HISAT2 not found in PATH" | tee -a HISAT2_pipeline_error.log; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] SAMtools not found in PATH" | tee -a HISAT2_pipeline_error.log; exit 1; }
command -v stringtie >/dev/null 2>&1 || { echo "[ERROR] StringTie not found in PATH" | tee -a HISAT2_pipeline_error.log; exit 1; }
command -v featureCounts >/dev/null 2>&1 || { echo "[ERROR] featureCounts not found in PATH" | tee -a HISAT2_pipeline_error.log; exit 1; }

# Set strand specificity for featureCounts:
# 0 = unstranded, 1 = stranded (forward), 2 = stranded (reverse)
# Most stranded RNA-seq protocols (e.g., Illumina TruSeq) use reverse-stranded libraries.
# Set strand specificity for featureCounts (2 = reverse-stranded)
STRAND=2


############################################################
# Dynamically locate HISAT2 helper scripts (Hybrid approach)
############################################################
echo "[INFO] Locating HISAT2 helper scripts" | tee -a HISAT2_pipeline_error.log

HISAT2_DIR=$(dirname "$(which hisat2)")
FALLBACK_DIR="/work/jdllab/Shared_Bioinformatics_Tools/hisat2_helpers"

# Function to resolve script path with fallback
resolve_script() {
  local script_name="$1"
  local dynamic_path
  dynamic_path=$(find "$HISAT2_DIR/.." -name "$script_name" 2>/dev/null | head -n 1)

  if [[ -f "$dynamic_path" ]]; then
    echo "$dynamic_path"
  elif [[ -f "$FALLBACK_DIR/$script_name" ]]; then
    echo "[WARNING] $script_name not found dynamically. Using fallback: $FALLBACK_DIR/$script_name" | tee -a HISAT2_pipeline_error.log
    echo "$FALLBACK_DIR/$script_name"
  else
    echo "[ERROR] Could not locate $script_name in either dynamic or fallback paths." >&2 | tee -a HISAT2_pipeline_error.log
    return 1
  fi
}

# Resolve each script
PY_HISAT2_SNP_SCRIPT=$(resolve_script "hisat2_extract_snps_haplotypes_VCF.py")
echo "[DEBUG] SNP script resolved to: $PY_HISAT2_SNP_SCRIPT" | tee -a HISAT2_pipeline_error.log

PY_HISAT2_SPLICE=$(resolve_script "hisat2_extract_splice_sites.py")
echo "[DEBUG] Splice site script resolved to: $PY_HISAT2_SPLICE" | tee -a HISAT2_pipeline_error.log

PY_HISAT2_EXONS=$(resolve_script "hisat2_extract_exons.py")
echo "[DEBUG] Exon script resolved to: $PY_HISAT2_EXONS" | tee -a HISAT2_pipeline_error.log

# Echo resolved paths
echo "[INFO] Using SNP script: $PY_HISAT2_SNP_SCRIPT" | tee -a HISAT2_pipeline_error.log
echo "[INFO] Using splice site script: $PY_HISAT2_SPLICE" | tee -a HISAT2_pipeline_error.log
echo "[INFO] Using exon script: $PY_HISAT2_EXONS" | tee -a HISAT2_pipeline_error.log



##########################################
#Preflight checks

#Early Error Detection
#Improved Debugging
##########################################

# Check BASE_DIR_PATH
echo "[INFO] BASE_DIR_PATH is: '$BASE_DIR_PATH'" | tee -a HISAT2_pipeline_error.log


# Check if trimmed data directory exists
if [[ ! -d "$RNA_SEQ_TRIMMED_FASTQ_PATH" ]]; then
    echo "[ERROR] Trimmed data directory not found: $RNA_SEQ_TRIMMED_FASTQ_PATH" | tee -a HISAT2_pipeline_error.log
    exit 1
else
    echo "[OK] Trimmed data directory found: $RNA_SEQ_TRIMMED_FASTQ_PATH"
fi

# Check if genome directory exists
if [[ ! -d "$GENOME_DIR_PATH" ]]; then
  echo "[ERROR] Genome directory not found: $GENOME_DIR_PATH" | tee -a HISAT2_pipeline_error.log
  exit 1
else
  echo "[OK] Genome directory found: $GENOME_DIR_PATH"
fi


#Check if genome assembly file exists and is not compressed
echo "[INFO] Checking for genome assembly FASTA file" | tee -a HISAT2_pipeline_error.log
if [[ -f "${GENOME_ASSEMBLY_PATH}.fna" ]]; then
    GENOME_FASTA="${GENOME_ASSEMBLY_PATH}.fna"
    echo "[OK] Genome assembly FASTA file found: $GENOME_FASTA"
    echo "[OK] Genome assembly FASTA file is plain text (not compressed)."
elif [[ -f "${GENOME_ASSEMBLY_PATH}.fa" ]]; then
    GENOME_FASTA="${GENOME_ASSEMBLY_PATH}.fa"
    echo "[OK] Genome FASTA found: $GENOME_FASTA"
    echo "[OK] Genome assembly FASTA file is plain text (not compressed)."
elif [[ -f "${GENOME_ASSEMBLY_PATH}.fa.gz" ]]; then
    echo "[WARNING] Compressed genome FASTA found: ${GENOME_ASSEMBLY_PATH}.fa.gz" | tee -a HISAT2_pipeline_error.log
    echo "[ERROR] Please decompress the .fa.gz file before running the pipeline." | tee -a HISAT2_pipeline_error.log
    exit 1
elif [[ -f "${GENOME_ASSEMBLY_PATH}.fna.gz" ]]; then
    echo "[WARNING] Compressed genome FASTA found: ${GENOME_ASSEMBLY_PATH}.fna.gz" | tee -a HISAT2_pipeline_error.log
    echo "[ERROR] Please decompress the .fna.gz file before running the pipeline." | tee -a HISAT2_pipeline_error.log
    exit 1
else
    echo "[ERROR] No genome FASTA file (.fa, .fna) found at expected location." | tee -a HISAT2_pipeline_error.log
    exit 1
fi


# Check if annotation file exists and is not compressed
if [[ ! -f "$ANNOTATION" ]]; then
    echo "[ERROR] Annotation file not found: $ANNOTATION" | tee -a HISAT2_pipeline_error.log
    exit 1
else
    echo "[OK] Annotation file found: $ANNOTATION"
fi

if [[ "$ANNOTATION" == *.gz ]]; then
    echo "[ERROR] Annotation file appears to be compressed. Please use an uncompressed .gtf file." | tee -a HISAT2_pipeline_error.log
    exit 1
else
    echo "[OK] Annotation file is plain text (not compressed)."
fi


echo "[INFO] Double-checking selected HISAT2 index path: $HISAT2_GENOME_INDEX_PATH"
MISSING=0
for i in {1..8}; do
    INDEX_FILE="${HISAT2_GENOME_INDEX_PATH}.${i}.ht2"
    if [[ -f "$INDEX_FILE" ]]; then
        echo "[OK] HISAT2 index file found: $INDEX_FILE"
    else
        echo "[ERROR] HISAT2 index file missing: $INDEX_FILE" | tee -a HISAT2_pipeline_error.log
        MISSING=1
    fi
done

if [[ $MISSING -eq 1 ]]; then
    echo "[FATAL] One or more HISAT2 index files are missing. Please verify the index path." | tee -a HISAT2_pipeline_error.log
    exit 1
fi


# Check if paired-end FASTQ files exist
if ! compgen -G "$RNA_SEQ_TRIMMED_FASTQ_PATH/*_1_paired.fq.gz" > /dev/null; then
    echo "[ERROR] No paired-end FASTQ files found in $RNA_SEQ_TRIMMED_FASTQ_PATH" | tee -a HISAT2_pipeline_error.log
    exit 1
else
    echo "[OK] Paired-end FASTQ files found in $RNA_SEQ_TRIMMED_FASTQ_PATH"
fi


#############################################
# Check and build HISAT2 index if not present
#############################################

echo "[INFO] Starting HISAT2 Index check" | tee -a HISAT2_pipeline_error.log

# Check if HISAT2 index files exist; build them if missing
if [[ -f "${HISAT2_GENOME_INDEX_PATH}.1.ht2" && -f "${HISAT2_GENOME_INDEX_PATH}.2.ht2" && \
      -f "${HISAT2_GENOME_INDEX_PATH}.3.ht2" && -f "${HISAT2_GENOME_INDEX_PATH}.4.ht2" && \
      -f "${HISAT2_GENOME_INDEX_PATH}.5.ht2" && -f "${HISAT2_GENOME_INDEX_PATH}.6.ht2" && \
      -f "${HISAT2_GENOME_INDEX_PATH}.7.ht2" && -f "${HISAT2_GENOME_INDEX_PATH}.8.ht2" ]]; then
    echo "[OK] Pre-built HISAT2 index files found. Skipping index build." | tee -a HISAT2_pipeline_error.log
else
    echo "[INFO] Pre-built HISAT2 index files not found. Checking for genome FASTA to build index..." | tee -a HISAT2_pipeline_error.log

    if [[ -f "${GENOME_ASSEMBLY_PATH}.fna" ]]; then
        GENOME_FASTA="${GENOME_ASSEMBLY_PATH}.fna"
    elif [[ -f "${GENOME_ASSEMBLY_PATH}.fa" ]]; then
        GENOME_FASTA="${GENOME_ASSEMBLY_PATH}.fa"
    else
        echo "[ERROR] Neither .fna nor .fa genome FASTA found. Cannot proceed." | tee -a HISAT2_pipeline_error.log
        exit 1
    fi
    echo "[OK] Genome FASTA found: $GENOME_FASTA" | tee -a HISAT2_pipeline_error.log

    mkdir -p "$HISAT2_GENOME_INDEX_PATH"
    echo "[INFO] Building HISAT2 index from genome FASTA..." | tee -a HISAT2_pipeline_error.log

    SNP_TXT="$HISAT2_GENOME_INDEX_PATH/snps.txt"
    HAPLOTYPE_TXT="$HISAT2_GENOME_INDEX_PATH/haplotype.txt"
    SNP_EXTRACTION_LOG="$HISAT2_GENOME_INDEX_PATH/snp_extraction.log"

    if [[ -f "$VCF_FILE_PATH" ]]; then
        echo "[INFO] SNP VCF file found: $VCF_FILE" | tee -a HISAT2_pipeline_error.log
        echo "[INFO] Building genome_snp_tran index with SNPs and transcript annotations..." | tee -a HISAT2_pipeline_error.log

        echo "[INFO] Extracting SNPs and haplotypes from VCF..." | tee -a "$SNP_EXTRACTION_LOG"
        "$PY_HISAT2_SNP_SCRIPT" "$GENOME_FASTA" "$VCF_FILE_PATH" "$HISAT2_GENOME_INDEX_PATH" 2>> "$SNP_EXTRACTION_LOG"
        echo "[INFO] SNP extraction and haplotype completed." | tee -a "$SNP_EXTRACTION_LOG"

        echo "[INFO] Building HISAT2 SNP- and transcript-aware index..." | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG"
        if ! hisat2-build \
            --snp "$SNP_TXT" \
            --haplotype "$HAPLOTYPE_TXT" \
            "$GENOME_FASTA" \
            "$HISAT2_GENOME_INDEX_PATH" \
            2>> "$HISAT2_GENOME_INDEX_BUILD_LOG"; then
            echo "[ERROR] genome_snp_tran index build failed. Falling back to genome_tran index..." | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG" | tee -a HISAT2_pipeline_error.log
            INDEX_SUBDIR="genome_tran"
            HISAT2_GENOME_INDEX_PATH="$GENOME_DIR_PATH/$INDEX_SUBDIR/$INDEX_SUBDIR"
            echo "[INFO] HISAT2_GENOME_INDEX_PATH reset to: $HISAT2_GENOME_INDEX_PATH" | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG"
            if ! hisat2-build "$GENOME_FASTA" "$HISAT2_GENOME_INDEX_PATH" 2>> "$HISAT2_GENOME_INDEX_BUILD_LOG"; then
                echo "[FATAL] Fallback genome_tran index build failed." | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG"
                exit 1
            fi
            echo "[INFO] HISAT2 genome_tran index build completed." | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG"
        fi
    else
        echo "[WARNING] SNP VCF file not found. Building genome_tran index without SNPs..." | tee -a HISAT2_pipeline_error.log
        if ! hisat2-build "$GENOME_FASTA" "$HISAT2_GENOME_INDEX_PATH" 2>> "$HISAT2_GENOME_INDEX_BUILD_LOG"; then
            echo "[FATAL] genome_tran index build failed." | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG"
            exit 1
        fi
        echo "[INFO] HISAT2 genome_tran index build completed." | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG"
    fi
fi


########################################
#RNAseq Mapping
########################################

#HISAT2: Here HISAT2 mapping is embedded inside a while loop that processes each
# paired-end read file.

#HISAT2 Usage: 
#  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]

#  <ht2-idx>  Index filename prefix (minus trailing .X.ht2).
#  <m1>       Files with #1 mates, paired with files in <m2>.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <m2>       Files with #2 mates, paired with files in <m1>.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <r>        Files with unpaired reads.
#             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#  <SRA accession number>        Comma-separated list of SRA accession numbers, e.g. --sra-acc SRR553653,SRR553654.
#  <sam>      File for SAM output (default: stdout)

#  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
#  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.:


#Comments organized to match lines below:
# -p # can use up to 8 threads to speed performance
# --dta option: Report alignments tailored for transcript assemblers including StringTie
# -x option: Path to the HISAT2 index for the reference genome (see above)
# -1 option: Path to the first mate (paired-end) reads
# -2 option: Path to the second mate (paired-end) reads
# -S option: Output file in SAM format

##############################################
echo "[INFO] Starting RNAseq Mapping Loop" | tee -a "$HISAT2_GENOME_INDEX_BUILD_LOG" 

# Extract splice sites and exons
echo "Extracting splice sites and exons from annotation..." | tee -a HISAT2_pipeline_error.log
echo "[DEBUG] Using annotation file: $ANNOTATION" | tee -a HISAT2_pipeline_error.log


# Extract splice sites and exons from annotation for HISAT2
$PY_HISAT2_SPLICE "$ANNOTATION" > "$RNA_SEQ_MAPPING_OUTPUT_DIR/splice_sites.txt" 2>> HISAT2_pipeline_error.log
$PY_HISAT2_EXONS "$ANNOTATION" > "$RNA_SEQ_MAPPING_OUTPUT_DIR/exons.txt" 2>> HISAT2_pipeline_error.log


# Create a list to store BAM files for featureCounts
BAM_LIST="$RNA_SEQ_MAPPING_OUTPUT_DIR/bam_files.txt"

# Ensure the output directory exists
mkdir -p "$(dirname "$BAM_LIST")"

# Clear the BAM list file and check for write errors
if ! : > "$BAM_LIST"; then
  echo "[ERROR] Could not write to $BAM_LIST" | tee -a HISAT2_pipeline_error.log
  exit 1
fi


############################################
#HISAT2 Loop
############################################
echo "[INFO] Starting HISAT2 StringTie Loop" | tee -a HISAT2_pipeline_error.log

# Loop through all paired-end FASTQ files for alignment and quantification
echo "[INFO] Looping through all paired-end trimmed files" | tee -a HISAT2_pipeline_error.log
echo "[DEBUG] Starting HISAT2 alignment loop..." | tee -a HISAT2_pipeline_error.log

# Validate and List Paired-End FASTQ Files
found_1=$(compgen -G "$RNA_SEQ_TRIMMED_FASTQ_PATH/*_1_paired.fq.gz")
found_2=$(compgen -G "$RNA_SEQ_TRIMMED_FASTQ_PATH/*_2_paired.fq.gz")

if [[ -n "$found_1" && -n "$found_2" ]]; then
  echo "[INFO] Found the following *_1_paired.fq.gz files:"
  printf "%s\n" $found_1
  echo "[INFO] Found the following *_2_paired.fq.gz files:"
  printf "%s\n" $found_2
else
  echo "[ERROR] Missing one or both paired-end FASTQ files in $RNA_SEQ_TRIMMED_FASTQ_PATH" | tee -a HISAT2_pipeline_error.log
  [[ -z "$found_1" ]] && echo "[ERROR] No *_1_paired.fq.gz files found." | tee -a HISAT2_pipeline_error.log
  [[ -z "$found_2" ]] && echo "[ERROR] No *_2_paired.fq.gz files found." | tee -a HISAT2_pipeline_error.log
  exit 1
fi


# === HISAT2 + StringTie Loop ===
for R1 in $found_1; do
    R2="${R1/_1_paired.fq.gz/_2_paired.fq.gz}"
    base=$(basename "$R1" _1_paired.fq.gz)

    SAM_OUT="$RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_hisat2.sam"
    BAM_OUT="$RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_hisat2.bam"
    SORTED_BAM="$RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_hisat2.sorted.bam"
    SUMMARY_OUT="$RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_hisat2_summary.txt"
    GTF_OUT="$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_stringtie.gtf"

    echo "[INFO] Processing $base..."

    if [[ ! -f "$SAM_OUT" ]]; then
        echo "[INFO] Mapping $base with HISAT2..."
        time hisat2 -p 8 --dta --summary-file "$SUMMARY_OUT" \
            --known-splicesite-infile "$RNA_SEQ_MAPPING_OUTPUT_DIR/splice_sites.txt" \
            -x "$HISAT2_GENOME_INDEX_PATH" -1 "$R1" -2 "$R2" -S "$SAM_OUT" 2>> HISAT2_pipeline_error.log
    else
        echo "[INFO] Skipping HISAT2 mapping for $base (already done)."
    fi

    if [[ ! -f "$BAM_OUT" ]]; then
        echo "[INFO] Converting SAM to BAM for $base..."
        samtools view -@ 4 -bS "$SAM_OUT" > "$BAM_OUT" 2>> HISAT2_pipeline_error.log
    else
        echo "[INFO] Skipping BAM conversion for $base (already done)."
    fi

    if [[ ! -f "$SORTED_BAM" ]]; then
        echo "[INFO] Sorting BAM for $base..."
        samtools sort -@ 4 "$BAM_OUT" -o "$SORTED_BAM" 2>> HISAT2_pipeline_error.log
    else
        echo "[INFO] Skipping BAM sorting for $base (already done)."
    fi

    if [[ ! -f "$SORTED_BAM.bai" ]]; then
        echo "[INFO] Indexing sorted BAM for $base..."
        samtools index "$SORTED_BAM" 2>> HISAT2_pipeline_error.log
    else
        echo "[INFO] Skipping BAM indexing for $base (already done)."
    fi

    echo "[INFO] Running StringTie for $base..."
    if [[ ! -f "$GTF_OUT" ]]; then
        stringtie "$SORTED_BAM" -p 8 -o "$GTF_OUT" -G "$ANNOTATION" --rf 2>> HISAT2_pipeline_error.log
        if [[ $? -ne 0 ]]; then
            echo "[WARNING] Reference-guided StringTie failed for $base. Retrying de novo..." | tee -a HISAT2_pipeline_error.log
            stringtie "$SORTED_BAM" -p 8 -o "$GTF_OUT" --rf 2>> HISAT2_pipeline_error.log
        fi
    else
        echo "[INFO] Skipping StringTie for $base (already done)."
    fi

    echo "$SORTED_BAM" >> "$RNA_SEQ_MAPPING_OUTPUT_DIR/bam_files.txt"
done

# Sort and deduplicate BAM list
sort -u "$BAM_LIST" -o "$BAM_LIST"


# Merge all GTFs into a unified transcriptome
MERGE_LIST="$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR/gtf_list.txt"
MERGED_GTF="$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR/stringtie_merged.gtf"

# Create list of all GTFs
find "$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR" -type f -name "*.gtf" > "$MERGE_LIST"

# Count how many GTFs were found
GTF_COUNT=$(wc -l < "$MERGE_LIST")
echo "[INFO] Found $GTF_COUNT GTF files for merging." | tee -a HISAT2_pipeline_error.log

# Check if the list is non-empty before attempting merge
if [[ "$GTF_COUNT" -eq 0 ]]; then
    echo "[ERROR] No GTF files found to merge. Skipping merge step." | tee -a HISAT2_pipeline_error.log
else
    # === Begin validation addition ===
    if [[ -f "$MERGED_GTF" && ! -s "$MERGED_GTF" ]]; then
        echo "[WARNING] Merged GTF exists but is empty. It will be regenerated." | tee -a HISAT2_pipeline_error.log
        rm "$MERGED_GTF"
    fi
    # === End validation addition ===

    if [[ ! -f "$MERGED_GTF" ]]; then
        echo "[INFO] $(date) - Starting GTF merge" | tee -a HISAT2_pipeline_error.log
        echo "Merging all GTFs into a unified transcriptome..."
        echo "[CMD] Running: stringtie --merge -p 8 -G $ANNOTATION -o $MERGED_GTF $MERGE_LIST"
        stringtie --merge -p 8 -G "$ANNOTATION" -o "$MERGED_GTF" "$MERGE_LIST" 2>> HISAT2_pipeline_error.log
        echo "[INFO] $(date) - Finished GTF merge with exit code $?" | tee -a HISAT2_pipeline_error.log
    else
        echo "Skipping GTF merge (already done)."
    fi
fi


# Re-quantify expression using the merged GTF
for R1 in $found_1; do
    base=$(basename "$R1" _1_paired.fq.gz)
    SORTED_BAM="$RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_hisat2.sorted.bam"
    BALLGOWN_DIR="$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR/${base}_ballgown"
    mkdir -p "$BALLGOWN_DIR"

    if [[ ! -f "$BALLGOWN_DIR/${base}.gtf" ]]; then
        echo "[INFO] $(date) - Starting re-quantification for $base" | tee -a HISAT2_pipeline_error.log
        echo "Re-quantifying $base using merged GTF..."
        echo "[CMD] Running: stringtie $SORTED_BAM -e -B -p 8 -G $MERGED_GTF -o $BALLGOWN_DIR/${base}.gtf"
        stringtie "$SORTED_BAM" -e -B -p 8 -G "$MERGED_GTF" -o "$BALLGOWN_DIR/${base}.gtf" 2>> HISAT2_pipeline_error.log
        echo "[INFO] $(date) - Finished re-quantification for $base with exit code $?" | tee -a HISAT2_pipeline_error.log
    else
        echo "Skipping re-quantification for $base (already done)."
    fi
done


##############################################
# Run featureCounts on all sorted BAM files
##############################################

# Run featureCounts on all sorted BAM files
if [[ ! -f "$COUNTS_OUT/gene_counts.txt" ]]; then
    echo "[INFO] $(date) - Running featureCounts on all BAM files..." | tee -a HISAT2_pipeline_error.log

    if [[ ! -s "$BAM_LIST" ]]; then
        echo "[ERROR] No BAM files listed in $BAM_LIST. Skipping featureCounts." | tee -a HISAT2_pipeline_error.log
        exit 1
    fi

    BAM_COUNT=$(wc -l < "$BAM_LIST")
    echo "[INFO] $BAM_COUNT BAM files listed for featureCounts." | tee -a HISAT2_pipeline_error.log

    featureCounts -T 8 \
        -a "$ANNOTATION" \
        -F GTF \
        -o "$COUNTS_OUT/gene_counts_featureCounts.txt" \
        -p -B -C \
        -t exon -g gene_id \
        -s $STRAND -O -M --primary -Q 10 \
        $(cat "$BAM_LIST") \
        2>> HISAT2_pipeline_error.log

    if [[ ! -s "$COUNTS_OUT/gene_counts.txt" ]]; then
        echo "[ERROR] featureCounts output file is empty. Please check for errors." | tee -a HISAT2_pipeline_error.log
        exit 1
    else
        echo "[INFO] featureCounts completed successfully and output is non-empty." | tee -a HISAT2_pipeline_error.log

#####################################################################
# === Run prepDE.py3.py to extract gene and transcript count matrices ===
#####################################################################

echo "[INFO] Running prepDE.py3.py to generate count matrices..."
PREPDE_SCRIPT="/work/jdllab/Shared_Bioinformatics_Tools/prepDE.py3"
BALLGOWN_DIR="$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR"
GENE_MATRIX="$COUNTS_OUT/gene_counts_prepDE.csv"
TX_MATRIX="$COUNTS_OUT/transcript_counts_prepDE.csv"
python3 "$PREPDE_SCRIPT" -i "$BALLGOWN_DIR" -g "$GENE_MATRIX" -t "$TX_MATRIX" 2>> HISAT2_pipeline_error.log
if [[ -s "$GENE_MATRIX" && -s "$TX_MATRIX" ]]; then
  echo "[INFO] prepDE.py3 completed successfully."
else
  echo "[ERROR] prepDE.py3 failed or output files are empty." | tee -a HISAT2_pipeline_error.log
  exit 1
fi
    fi
else
    echo "[INFO] $(date) - Skipping featureCounts (already done)." | tee -a HISAT2_pipeline_error.log
fi

###############################
# === Final Output Summary ===
###############################
echo "Final output summary:"
echo "[SUMMARY] Pipeline completed at: $(date)"
echo "[SUMMARY] featureCounts gene matrix: $COUNTS_OUT/gene_counts_featureCounts.txt"
echo "[SUMMARY] prepDE gene matrix: $GENE_MATRIX"
echo "[SUMMARY] prepDE transcript matrix: $TX_MATRIX"
echo "[SUMMARY] Total samples processed: $(wc -l < "$BAM_LIST")"
echo "[SUMMARY] Total GTFs merged: $GTF_COUNT"
echo "Merged GTF: $MERGED_GTF"
echo "Ballgown directories:"
find "$GTF_RNA_SEQ_MAPPING_OUTPUT_DIR" -type d -name "*_ballgown"

echo -e "\n** Script ended on $(date) **\n"
echo "Done!"

############################################################
# Annotation: Understanding Count Matrices in This Pipeline
#
# This pipeline generates two types of gene expression count matrices:
#
# 1. featureCounts Output
#    - File: gene_counts.txt
#    - Source: Counts reads overlapping exons in BAM files.
#    - Tool: featureCounts (Subread)
#    - Level: Gene-level only
#    - Use when:
#        * You want raw gene counts for DESeq2, edgeR, or limma-voom.
#        * You prefer fast, robust, strand-specific quantification.
#        * You are working with well-annotated genomes.
#
# 2. prepDE.py Output
#    - Files: gene_count_matrix_prepDE.csv, transcript_count_matrix_prepDE.csv
#    - Source: Extracted from StringTie quantification (Ballgown directories)
#    - Tool: prepDE.py
#    - Level: Gene- and transcript-level
#    - Use when:
#        * You need transcript-level or isoform-specific DE.
#        * You are using StringTie’s reference-guided or de novo assembly.
#
# Which Should You Use for DESeq2?
# | Goal                          | Recommended Matrix                  |
# |-------------------------------|-------------------------------------|
# | Standard gene-level DE        | gene_counts.txt (featureCounts)     |
# | Transcript-level DE           | transcript_count_matrix_prepDE.csv  |
# | Gene-level DE from StringTie  | gene_count_matrix_prepDE.csv        |
# | Isoform-specific analysis     | transcript_count_matrix_prepDE.csv  |
#
# ⚠️ Do not mix featureCounts and prepDE.py outputs in the same DESeq2 analysis.
############################################################
