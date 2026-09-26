#!/usr/bin/env bash
#
# pipeline_common.sh -- shared helpers for the PAX6 RNA-seq pipeline.
#
# James D. Lauderdale, PhD; Department of Cellular Biology,
# University of Georgia, Athens, GA 30602, USA
# Study: "Nerve remodeling in a Pax6 model of keratopathy"
#
# Sourced by every pipeline script. Contains no pipeline logic of its own:
# paths, study/batch definitions and tool versions live in data files
# (studies.tsv, module_versions.txt) so that adding data never means editing
# a script. Not executable on its own.

set -euo pipefail

# SLURM stages the batch script to a spool directory on the compute node, so
# dirname "$0" does NOT point at the submit directory. Prefer SLURM_SUBMIT_DIR.
PIPELINE_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
MANIFEST="${PAX6_MANIFEST:-$PIPELINE_DIR/studies.tsv}"
VERSION_FILE="${PAX6_VERSION_FILE:-$PIPELINE_DIR/module_versions.txt}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

PROJECT_ROOT="${PAX6_ROOT:-/work/jdllab/PAX6_RNAseq}"
ASSEMBLY_FILE="${PAX6_ASSEMBLY_FILE:-$PIPELINE_DIR/assemblies.tsv}"


log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die()  { printf '[FATAL] %s\n' "$*" >&2; exit 1; }
warn() { printf '[WARN]  %s\n' "$*" >&2; }

# ---------------------------------------------------------------------------
# Module loading, pinned
# ---------------------------------------------------------------------------
# The earlier scripts either hard-coded a module or loaded "the latest
# available version" at run time, and neither recorded which version ran.
# Subread changed the meaning of -p at 2.0.2; the old counting command was
# written for the earlier meaning and ran, in both batches, under 2.0.6. It
# counted reads instead of fragments without a single error. Pinning, plus
# report_versions() in every log, makes a change of that kind visible. A pinned
# version that fails loudly beats a floating version that succeeds quietly.
#
# To move a pin: change it deliberately, then RE-COUNT EVERY BATCH in the
# affected study. Never edit the pin merely to make a job start.
load_pinned() {
  local want name spec
  [[ -f "$VERSION_FILE" ]] || die "Version file not found: $VERSION_FILE"
  module purge 2>/dev/null || true
  for want in "$@"; do
    spec=$(awk -v n="$want" '$1 == n && $0 !~ /^#/ { print $2; exit }' "$VERSION_FILE")
    [[ -n "$spec" ]] || die "No pinned version for '$want' in $VERSION_FILE"
    name="$want/$spec"
    if ! ml "$name" 2>/dev/null; then
      echo "[FATAL] Module not available: $name" >&2
      echo "        Versions present on this system:" >&2
      module -t avail "$want" 2>&1 | sed 's/^/          /' >&2
      echo "        Either keep that version, or change the pin and re-count." >&2
      exit 1
    fi
    log "loaded $name"
  done
}

report_versions() {
  local tool
  for tool in fastqc multiqc hisat2 samtools featureCounts; do
    command -v "$tool" >/dev/null 2>&1 || continue
    case "$tool" in
      hisat2)        printf '  %-14s %s\n' hisat2        "$(hisat2 --version 2>&1 | head -n1)" ;;
      samtools)      printf '  %-14s %s\n' samtools      "$(samtools --version 2>&1 | head -n1)" ;;
      featureCounts) printf '  %-14s %s\n' featureCounts "$(featureCounts -v 2>&1 | grep -o 'v[0-9.]*' | head -n1)" ;;
      fastqc)        printf '  %-14s %s\n' fastqc        "$(fastqc --version 2>&1 | head -n1)" ;;
      multiqc)       printf '  %-14s %s\n' multiqc       "$(multiqc --version 2>&1 | head -n1)" ;;
    esac
  done
}

# ---------------------------------------------------------------------------
# Manifest: study / batch / layout / strand / trimmed_dir / mapping_dir
# ---------------------------------------------------------------------------
_mf() { awk -F'\t' '$0 !~ /^#/ && NF >= 6' "$MANIFEST"; }

manifest_studies()  { _mf | awk -F'\t' '{ print $1 }' | sort -u; }
manifest_batches()  { _mf | awk -F'\t' -v s="$1" '$1 == s { print $2 }'; }
manifest_all_keys() { _mf | awk -F'\t' '{ print $1"/"$2 }'; }

# manifest_field <study> <batch> <col>
manifest_field() {
  _mf | awk -F'\t' -v s="$1" -v b="$2" -v c="$3" \
    '$1 == s && $2 == b { gsub(/^ +| +$/, "", $c); print $c; exit }'
}

require_manifest() {
  [[ -f "$MANIFEST" ]] || die "Manifest not found: $MANIFEST"
  local n dupes
  n=$(_mf | wc -l)
  (( n > 0 )) || die "No rows in $MANIFEST"
  dupes=$(manifest_all_keys | sort | uniq -d)
  [[ -z "$dupes" ]] || die "Duplicate study/batch in $MANIFEST: $dupes"
}

require_study() {
  manifest_studies | grep -qx "$1" || {
    echo "[FATAL] Study '$1' is not in $MANIFEST. Studies defined:" >&2
    manifest_studies | sed 's/^/          /' >&2
    exit 1
  }
}


# ---------------------------------------------------------------------------
# Assembly selection
# ---------------------------------------------------------------------------
# Assembly is a property of a COUNT MATRIX, not of a batch. BAMs go to
# <mapping_root>/<assembly>/, so the same libraries can exist against several
# assemblies at once; counting reads one assembly's directory and records
# which, so a matrix cannot contain a mixture.
_af() { awk -F'\t' '$0 !~ /^#/ && NF >= 6' "$ASSEMBLY_FILE"; }

assembly_list()  { _af | awk -F'\t' '{ print $1 }'; }
assembly_field() { _af | awk -F'\t' -v a="$1" -v c="$2" '$1 == a { gsub(/^ +| +$/, "", $c); print $c; exit }'; }

resolve_assembly() {
  [[ -f "$ASSEMBLY_FILE" ]] || die "Assembly file not found: $ASSEMBLY_FILE"
  ASSEMBLY="${PAX6_ASSEMBLY:-$(_af | head -n1 | cut -f1)}"
  [[ -n "$ASSEMBLY" ]] || die "No assemblies defined in $ASSEMBLY_FILE"
  if ! assembly_list | grep -qx "$ASSEMBLY"; then
    echo "[FATAL] Unknown assembly '$ASSEMBLY'. Defined in $ASSEMBLY_FILE:" >&2
    assembly_list | sed 's/^/          /' >&2
    exit 1
  fi
  local ref
  ref=$(assembly_field "$ASSEMBLY" 2)
  GENOME_DIR="${PAX6_GENOME_DIR:-$ref}"
  GENOME_FASTA="$GENOME_DIR/$(assembly_field "$ASSEMBLY" 3)"
  ANNOTATION="${PAX6_ANNOTATION:-$GENOME_DIR/$(assembly_field "$ASSEMBLY" 4)}"
  HISAT2_INDEX="${PAX6_HISAT2_INDEX:-$GENOME_DIR/$(assembly_field "$ASSEMBLY" 5)}"
  INDEX_TYPE=$(assembly_field "$ASSEMBLY" 6)
}
resolve_assembly

# Where this batch's BAMs live FOR THE SELECTED ASSEMBLY.
mapping_dir_for() { printf '%s/%s\n' "$(manifest_field "$1" "$2" 6)" "$ASSEMBLY"; }

# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------
# One convention: <sample_id>_hisat2.sorted.bam. Study and batch identity live
# in the manifest and in sample_map.tsv, never in a filename. The trailing
# strips handle BAMs already on disk under the older 2026 convention.
sample_id_from_bam() {
  local b; b=$(basename "$1")
  b="${b%_hisat2.sorted.bam}"; b="${b%.sorted.bam}"; b="${b%.bam}"; b="${b%.sam}"
  b="${b%_Novogene_*}"; b="${b%_hisat2}"
  printf '%s\n' "$b"
}

sample_from_r1() {
  local b; b=$(basename "$1")
  b="${b%_1_paired.fq.gz}"; b="${b%_1.fastq.gz}"; b="${b%_1.fq.gz}"
  b="${b%_trimmed.fq.gz}";  b="${b%.fastq.gz}";   b="${b%.fq.gz}"
  printf '%s\n' "$b"
}

# -L throughout: directories and files in this project are frequently symlinks
# (the August 2026 reorganisation linked legacy data into the tree). Without -L, `find`
# tests the link itself, `-type f` is false, and the search silently returns
# nothing -- which looks exactly like "there is no data here".
# Optional second argument: search depth (default 1). Vendor deliveries often
# nest one directory per sample, so qc_trim_batch.sh searches deeper.
find_r1() {
  find -L "$1" -maxdepth "${2:-1}" -type f \
    \( -name "*_1_paired.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort
}

# Single-end FASTQs: anything that is not a member of a pair.
find_se() {
  find -L "$1" -maxdepth "${2:-1}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
    ! -name "*_1.fastq.gz" ! -name "*_2.fastq.gz" \
    ! -name "*_1.fq.gz"    ! -name "*_2.fq.gz" \
    ! -name "*_paired.fq.gz" ! -name "*_unpaired.fq.gz" | sort
}

mate_of() {
  local r1="$1" r2
  case "$r1" in
    *_1_paired.fq.gz) r2="${r1/_1_paired.fq.gz/_2_paired.fq.gz}" ;;
    *_1.fastq.gz)     r2="${r1/_1.fastq.gz/_2.fastq.gz}" ;;
    *_1.fq.gz)        r2="${r1/_1.fq.gz/_2.fq.gz}" ;;
    *) return 1 ;;
  esac
  [[ -f "$r2" ]] || return 1
  printf '%s\n' "$r2"
}

# A BAM is usable only if complete AND indexed. quickcheck catches the
# truncated file a killed job leaves behind, otherwise invisible until
# featureCounts reports nonsense.
bam_ok() {
  local bam="$1"
  [[ -s "$bam" ]] || return 1
  [[ -s "$bam.bai" || -s "${bam%.bam}.bai" ]] || return 1
  samtools quickcheck "$bam" 2>/dev/null || return 1
  return 0
}

# Measured, not assumed: PE if any alignment carries the paired flag (0x1).
detect_layout() {
  local n
  n=$(samtools view -c -f 1 "$1" 2>/dev/null) || return 1
  [[ "${n:-0}" -gt 0 ]] && printf 'PE\n' || printf 'SE\n'
}
