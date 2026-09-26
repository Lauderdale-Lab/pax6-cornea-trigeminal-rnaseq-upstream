# shellcheck shell=bash
#
# lib/common.sh -- shared functions for the PAX6 RNA-seq pipeline.
#
# Sourced by bin/pax6 and by every job in slurm/. Holds no analysis logic of
# its own. Paths, datasets, genome builds, tool versions and parameters all
# live in config/, so adding data never means editing code.

set -euo pipefail

# PAX6_HOME is the repository root. bin/pax6 exports it to every job, because
# SLURM copies a batch script to a spool directory and its own path is useless.
if [[ -z "${PAX6_HOME:-}" ]]; then
  PAX6_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export PAX6_HOME
CONFIG_DIR="$PAX6_HOME/config"

# shellcheck source=../config/site.env
source "$CONFIG_DIR/site.env"
export PAX6_ROOT PAX6_RUNS PAX6_KEEP
DATA_DIR="$PAX6_ROOT/data"
REF_DIR="$PAX6_ROOT/reference"
RUNS_DIR="$PAX6_RUNS"
KEEP_DIR="$PAX6_KEEP"

export THREADS="${SLURM_CPUS_PER_TASK:-4}"

# ---------------------------------------------------------------------------
# Messages
# ---------------------------------------------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
warn() { printf '[WARN]  %s\n' "$*" >&2; }
die()  { printf '[FATAL] %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Tab-separated config tables, read by COLUMN NAME
# ---------------------------------------------------------------------------
# Comment lines (#) and blank lines are ignored. The first remaining line is
# the header. Columns are looked up by name, so adding or reordering columns
# never breaks a caller.

# tsv_rows <file> : data rows only, header removed
tsv_rows() {
  awk -F'\t' '/^[[:space:]]*#/ || /^[[:space:]]*$/ { next } !seen++ { next } { print }' "$1"
}

# tsv_col <file> <column> : 1-based index of a named column, or die
tsv_col() {
  local idx
  idx=$(awk -F'\t' -v c="$2" '/^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
          { for (i = 1; i <= NF; i++) if ($i == c) { print i; exit } exit }' "$1")
  [[ -n "$idx" ]] || die "Column '$2' not found in $1"
  printf '%s\n' "$idx"
}

# tsv_get <file> <key_value> <column> : value of <column> in the row whose
# FIRST column equals <key_value>; empty if no such row
tsv_get() {
  local c; c=$(tsv_col "$1" "$3")
  tsv_rows "$1" | awk -F'\t' -v k="$2" -v c="$c" '$1 == k { print $c; exit }'
}

# ---------------------------------------------------------------------------
# Datasets and studies
# ---------------------------------------------------------------------------
DATASETS_TSV="$CONFIG_DIR/datasets.tsv"

dataset_list()   { tsv_rows "$DATASETS_TSV" | cut -f1; }
dataset_field()  { tsv_get "$DATASETS_TSV" "$1" "$2"; }
study_list()     { local c; c=$(tsv_col "$DATASETS_TSV" study); tsv_rows "$DATASETS_TSV" | cut -f"$c" | sort -u; }
study_datasets() { local c; c=$(tsv_col "$DATASETS_TSV" study); tsv_rows "$DATASETS_TSV" | awk -F'\t' -v s="$1" -v c="$c" '$c == s { print $1 }'; }

require_dataset() {
  dataset_list | grep -qxF "$1" || die "Dataset '$1' is not in $DATASETS_TSV. Defined: $(dataset_list | tr '\n' ' ')"
}
require_study() {
  study_list | grep -qxF "$1" || die "Study '$1' is not in $DATASETS_TSV. Defined: $(study_list | tr '\n' ' ')"
}
raw_dir_for() { printf '%s/%s/raw\n' "$DATA_DIR" "$1"; }

# ---------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------
# A run is one complete processing of data with ONE pipeline commit, ONE
# genome build, ONE set of module versions and ONE set of parameters, all
# frozen into runs/<run>/ when the run is created. Outputs never mix runs.

run_dir() { printf '%s/%s\n' "$RUNS_DIR" "$1"; }

run_info() {  # run_info <run> <key>
  awk -F'\t' -v k="$2" '$1 == k { print $2; exit }' "$(run_dir "$1")/RUN_INFO.tsv"
}

# load_run <run> : check the run exists, then load its frozen settings
load_run() {
  local r="$1" d
  [[ -n "$r" ]] || die "No run given."
  d=$(run_dir "$r")
  [[ -f "$d/RUN_INFO.tsv" ]] || die "Run '$r' does not exist ($d). Create it with: bin/pax6 new-run $r"
  RUN="$r"; RUN_DIR="$d"
  MODULES_TSV="$RUN_DIR/modules.tsv"
  # shellcheck source=../config/parameters.env
  source "$RUN_DIR/parameters.env"
  resolve_assembly "$(run_info "$r" assembly)"
  export RUN RUN_DIR
}

# Where each kind of output lives inside a run
ds_dir()      { printf '%s/%s\n' "$RUN_DIR" "$1"; }
samples_tsv() { printf '%s/%s/samples.tsv\n' "$RUN_DIR" "$1"; }
trim_dir()    { printf '%s/%s/trimmed\n' "$RUN_DIR" "$1"; }
qc_dir()      { printf '%s/%s/qc\n' "$RUN_DIR" "$1"; }
align_dir()   { printf '%s/%s/align\n' "$RUN_DIR" "$1"; }
strand_dir()  { printf '%s/%s/strand\n' "$RUN_DIR" "$1"; }
counts_dir()  { printf '%s/counts/%s\n' "$RUN_DIR" "$1"; }

# ---------------------------------------------------------------------------
# Genome builds
# ---------------------------------------------------------------------------
ASSEMBLIES_TSV="$CONFIG_DIR/assemblies.tsv"

resolve_assembly() {
  ASSEMBLY="$1"
  tsv_rows "$ASSEMBLIES_TSV" | cut -f1 | grep -qxF "$ASSEMBLY" ||
    die "Unknown assembly '$ASSEMBLY'. Defined: $(tsv_rows "$ASSEMBLIES_TSV" | cut -f1 | tr '\n' ' ')"
  GENOME_DIR="$REF_DIR/$ASSEMBLY"
  GENOME_FASTA="$GENOME_DIR/$(tsv_get "$ASSEMBLIES_TSV" "$ASSEMBLY" genome_fasta)"
  ANNOTATION="$GENOME_DIR/$(tsv_get "$ASSEMBLIES_TSV" "$ASSEMBLY" annotation)"
  HISAT2_INDEX="$GENOME_DIR/$(tsv_get "$ASSEMBLIES_TSV" "$ASSEMBLY" index_prefix)"
  INDEX_TYPE="$(tsv_get "$ASSEMBLIES_TSV" "$ASSEMBLY" index_type)"
  export ASSEMBLY GENOME_DIR GENOME_FASTA ANNOTATION HISAT2_INDEX INDEX_TYPE
}

# ---------------------------------------------------------------------------
# Modules, pinned
# ---------------------------------------------------------------------------
# A pinned module that fails loudly beats a floating one that succeeds
# quietly: the counting defect behind the superseded matrices came from a
# command written for one Subread behaviour and run under another.
#
# load_modules <tool>... : purge, then load each tool's pinned module
load_modules() {
  local tsv="${MODULES_TSV:-$CONFIG_DIR/modules.tsv}" tool mod f
  if ! type module >/dev/null 2>&1; then
    for f in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh; do
      # shellcheck disable=SC1090
      [[ -f "$f" ]] && { source "$f"; break; }
    done
  fi
  type module >/dev/null 2>&1 || die "The 'module' command is not available in this shell."
  module purge >/dev/null 2>&1 || true
  for tool in "$@"; do
    mod=$(tsv_get "$tsv" "$tool" module)
    [[ -n "$mod" ]] || die "No pinned module for '$tool' in $tsv"
    if ! module load "$mod" >/dev/null 2>&1; then
      echo "[FATAL] Module not available: $mod" >&2
      module -t avail "${mod%%/*}" 2>&1 | sed 's/^/        available: /' >&2 || true
      exit 1
    fi
  done
  log "modules: $*  ->  $(module -t list 2>&1 | grep -E "^($(IFS='|'; echo "$*"))/" | tr '\n' ' ')"
}

# ---------------------------------------------------------------------------
# FASTQ discovery
# ---------------------------------------------------------------------------
# Sample IDs come from file names: <sample>_1.fq.gz / <sample>_1.fastq.gz and
# their _2 mates, or <sample>.fq.gz for single-end. Vendors often nest one
# folder per sample, so the search goes several levels deep. find -L follows
# symlinks: without it a directory of links looks empty.

sample_from_fastq() {
  local b; b=$(basename "$1")
  b="${b%.gz}"; b="${b%.fastq}"; b="${b%.fq}"; b="${b%_[12]}"
  printf '%s\n' "$b"
}

mate_of() {
  local r1="$1" r2 dir base
  dir=$(dirname "$r1"); base=$(basename "$r1")
  r2="$dir/${base/_1.f/_2.f}"
  [[ "$r2" != "$r1" && -f "$r2" ]] || return 1
  printf '%s\n' "$r2"
}

# discover_samples <raw_dir> : prints "sample_id<TAB>r1<TAB>r2" (r2 is "-" for
# single-end), sorted by sample. Dies on a missing mate, a duplicate sample
# name, or a mixture of paired and single-end files.
discover_samples() {
  local raw="$1" depth="${PAX6_RAW_DEPTH:-4}" f s r2 n_pe n_se
  [[ -d "$raw" ]] || die "Raw directory not found: $raw"
  local pe se
  pe=$(find -L "$raw" -maxdepth "$depth" -type f \( -name '*_1.fq.gz' -o -name '*_1.fastq.gz' \) | sort)
  se=$(find -L "$raw" -maxdepth "$depth" -type f \( -name '*.fq.gz' -o -name '*.fastq.gz' \) \
         ! -name '*_1.f*q.gz' ! -name '*_2.f*q.gz' | sort)
  n_pe=$(grep -c . <<<"$pe" || true); n_se=$(grep -c . <<<"$se" || true)
  (( n_pe > 0 && n_se > 0 )) && die "$raw mixes paired-end (_1/_2) and single-end FASTQ files."
  (( n_pe + n_se > 0 )) || die "No FASTQ files (*.fq.gz, *.fastq.gz) under $raw"
  {
    if (( n_pe > 0 )); then
      while IFS= read -r f; do
        s=$(sample_from_fastq "$f")
        r2=$(mate_of "$f") || die "No _2 mate for $f"
        printf '%s\t%s\t%s\n' "$s" "$f" "$r2"
      done <<<"$pe"
    else
      while IFS= read -r f; do printf '%s\t%s\t-\n' "$(sample_from_fastq "$f")" "$f"; done <<<"$se"
    fi
  } | sort -k1,1 > "${TMPDIR:-/tmp}/pax6_discover.$$"
  local dups; dups=$(cut -f1 "${TMPDIR:-/tmp}/pax6_discover.$$" | uniq -d)
  if [[ -n "$dups" ]]; then
    rm -f "${TMPDIR:-/tmp}/pax6_discover.$$"
    die "Sample name(s) found more than once under $raw: $(tr '\n' ' ' <<<"$dups")"
  fi
  cat "${TMPDIR:-/tmp}/pax6_discover.$$"; rm -f "${TMPDIR:-/tmp}/pax6_discover.$$"
}

# sample_line <samples.tsv> <n> : the n-th sample (1-based), for array tasks
sample_line() {
  local line; line=$(sed -n "${2}p" "$1")
  [[ -n "$line" ]] || die "No sample number $2 in $1"
  printf '%s\n' "$line"
}

# ---------------------------------------------------------------------------
# Output verification
# ---------------------------------------------------------------------------
# A step is complete only when its product verifies, never because a file
# name exists: a job killed mid-write leaves a file that looks finished.

bam_ok() {
  [[ -s "$1" && -s "$1.bai" ]] || return 1
  samtools quickcheck "$1" 2>/dev/null
}

# gz_ok <file> : non-empty and a complete gzip stream
gz_ok() { [[ -s "$1" ]] && gzip -t "$1" 2>/dev/null; }

# ---------------------------------------------------------------------------
# featureCounts summaries
# ---------------------------------------------------------------------------
# assigned_fraction <summary> : Assigned / total, first sample column
assigned_fraction() {
  awk -F'\t' 'NR > 1 { t += $2; v[$1] = $2 } END { printf "%.4f\n", (t ? v["Assigned"] / t : 0) }' "$1"
}

# strand_call <a0> <a1> <a2> : strandedness from the assigned fractions at
# -s 0, 1 and 2. Prints 0, 1, 2 or AMBIGUOUS.
#   stranded:   the matching setting keeps most of -s 0, the other loses most
#   unstranded: -s 1 and -s 2 each keep roughly half of -s 0
strand_call() {
  awk -v a0="$1" -v a1="$2" -v a2="$3" 'BEGIN {
    if (a0 <= 0)                                           { print "AMBIGUOUS"; exit }
    r1 = a1 / a0; r2 = a2 / a0
    if (r2 > 0.6 && r1 < 0.25)                             print "2"
    else if (r1 > 0.6 && r2 < 0.25)                        print "1"
    else if (r1 > 0.3 && r1 < 0.7 && r2 > 0.3 && r2 < 0.7) print "0"
    else                                                   print "AMBIGUOUS"
  }'
}

# assignment_table <featureCounts summary> : one row per library with
# Assigned %, NoFeatures % and Assigned/(Assigned+NoFeatures) %. Column
# headers are reduced to sample IDs.
assignment_table() {
  awk -F'\t' '
    NR == 1 { for (i = 2; i <= NF; i++) { n = split($i, p, "/"); s = p[n]; sub(/(_hisat2\.sorted)?\.bam$/, "", s); name[i] = s } ; nf = NF; next }
    { for (i = 2; i <= NF; i++) { tot[i] += $i; v[$1, i] = $i } }
    END {
      print "sample_id\tassigned_pct\tnofeatures_pct\tassigned_ratio_pct"
      for (i = 2; i <= nf; i++) {
        a = v["Assigned", i]; f = v["Unassigned_NoFeatures", i]
        printf "%s\t%.1f\t%.1f\t%.1f\n", name[i], 100*a/tot[i], 100*f/tot[i], ((a+f) ? 100*a/(a+f) : 0)
      }
    }' "$1"
}

# print_table [file] : align a tab-separated table for reading (stdin if no file)
print_table() {
  awk -F'\t' '{ for (i = 1; i <= NF; i++) { c[NR, i] = $i; if (length($i) > w[i]) w[i] = length($i) } ; if (NF > nf) nf = NF }
    END { for (r = 1; r <= NR; r++) { line = ""; for (i = 1; i <= nf; i++) line = line sprintf("%-" w[i] + 2 "s", c[r, i]); sub(/ +$/, "", line); print line } }' "${1:-/dev/stdin}"
}

# ---------------------------------------------------------------------------
# Software versions, as the tools themselves report them
# ---------------------------------------------------------------------------
# Module names say what was asked for; this records what actually ran. Each
# job writes one small file; count.sbatch gathers them into the matrix's
# provenance and flags any tool that ran at more than one version.

tool_version() {
  local jar
  case "$1" in
    FastQC)      fastqc --version 2>&1 | head -n1 ;;
    Trimmomatic) jar=$(find "${EBROOTTRIMMOMATIC:-/nonexistent}" -maxdepth 1 -name 'trimmomatic*.jar' 2>/dev/null | head -n1)
                 [[ -n "$jar" ]] && printf 'Trimmomatic %s\n' "$(java -jar "$jar" -version 2>&1 | tail -n1)" ;;
    MultiQC)     multiqc --version 2>&1 | head -n1 ;;
    HISAT2)      hisat2 --version 2>&1 | head -n1 | sed 's|^.*/||' ;;
    SAMtools)    samtools --version 2>&1 | head -n1 ;;
    Subread)     printf 'featureCounts %s\n' "$(featureCounts -v 2>&1 | grep -o 'v[0-9][0-9.]*' | head -n1)" ;;
    *)           echo unknown ;;
  esac
}

# record_versions <step> <tool>... : after load_modules, inside a run
record_versions() {
  local step="$1" t d="$RUN_DIR/provenance/versions"; shift
  mkdir -p "$d"
  for t in "$@"; do
    printf '%s\t%s\t%s\t%s\t%s\n' "$t" "$(tsv_get "${MODULES_TSV:-$CONFIG_DIR/modules.tsv}" "$t" module)" \
      "$(tool_version "$t")" "$step" "$(hostname 2>/dev/null || echo unknown)"
  done > "$d/${step}.${PAX6_DATASET:-${PAX6_STUDY:-all}}.${SLURM_JOB_ID:-local}.${SLURM_ARRAY_TASK_ID:-0}.tsv"
}

# stats <file> <column> : "min max median" of a numeric column (header skipped)
stats() {
  local c; c=$(tsv_col "$1" "$2")
  tail -n +2 "$1" | cut -f"$c" | grep -E '^[0-9.]+$' | sort -g |
    awk '{ v[NR] = $1 } END { if (!NR) { print "NA NA NA"; exit }
          m = (NR % 2) ? v[(NR + 1) / 2] : (v[NR / 2] + v[NR / 2 + 1]) / 2
          printf "%s %s %s\n", v[1], v[NR], m }'
}

# ---------------------------------------------------------------------------
# Keeping results off scratch
# ---------------------------------------------------------------------------
# keep_results <run> [--with-bams] : copy a run's durable results from
# $PAX6_RUNS (scratch) to $PAX6_KEEP (/work). Always copied: run records,
# count matrices and provenance, per-dataset samples, checksums, strand
# results, QC reports and logs. BAMs only with --with-bams. Trimmed reads are
# never copied: they are regenerated from raw/ in about an hour per library.
keep_results() {
  local r="$1" with_bams="${2:-}" src dst f
  src=$(run_dir "$r"); dst="$KEEP_DIR/$r"
  [[ -d "$src" ]] || die "Run $r not found at $src"
  mkdir -p "$dst"
  # Nothing to do when runs are already kept where they are processed
  [[ "$(cd "$src" && pwd -P)" == "$(cd "$dst" && pwd -P)" ]] && { printf '%s\n' "$dst"; return 0; }
  ( cd "$src" && find . -type f ! -path '*/trimmed/*' ) | while IFS= read -r f; do
    if [[ "$with_bams" != --with-bams ]] && [[ "$f" == *.bam || "$f" == *.bam.bai ]]; then continue; fi
    # copy only what is new or changed, so repeating is cheap
    if [[ ! -e "$dst/$f" || "$src/$f" -nt "$dst/$f" ]]; then
      mkdir -p "$dst/$(dirname "$f")"; cp -p "$src/$f" "$dst/$f"
    fi
  done
  printf '%s\n' "$dst"
}
