#!/usr/bin/env bash
#
# slurm/multiqc.sbatch -- one MultiQC report per dataset, covering raw and
# trimmed FastQC, Trimmomatic and HISAT2. Runs after alignment finishes.
#
#SBATCH --job-name=pax6_multiqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
load_run "$PAX6_RUN"
DATASET="$PAX6_DATASET"
load_modules MultiQC
record_versions multiqc MultiQC
# -d -dd 1 prefixes each sample with its folder (raw | trimmed | align).
# Raw and trimmed FastQC reports otherwise share sample names and MultiQC
# would keep only one of each pair.
multiqc -f -d -dd 1 -o "$(ds_dir "$DATASET")" -n "multiqc_${DATASET}.html" \
  "$(qc_dir "$DATASET")" "$(align_dir "$DATASET")"
