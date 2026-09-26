#!/usr/bin/env bash
#
# tools/make_test_dataset.sh <source_dataset> <read_pairs> <sample>...
#
# Makes a small dataset for testing the whole pipeline in minutes: the first
# <read_pairs> reads of each named library, written to
# data/TEST_<source_dataset>/raw/. Then add the row it prints to
# config/datasets.tsv and run it like any other dataset.
#
#   tools/make_test_dataset.sh Lauderdale_2026 200000 A_CW3 A_CU3
#
# The first reads of a file are not a random sample, which is fine for
# checking that every step runs; it is not a basis for any biological number.

set -euo pipefail
PAX6_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
SRC="${1:?usage: make_test_dataset.sh <source_dataset> <read_pairs> <sample>...}"
N="${2:?read pairs}"; shift 2
(( $# > 0 )) || die "Name at least one sample."
require_dataset "$SRC"
OUT="$DATA_DIR/TEST_$SRC/raw"
mkdir -p "$OUT"
ALL=$(discover_samples "$(raw_dir_for "$SRC")")
for S in "$@"; do
  line=$(awk -F'\t' -v s="$S" '$1 == s' <<<"$ALL")
  [[ -n "$line" ]] || die "Sample $S not found in $SRC"
  IFS=$'\t' read -r _ R1 R2 <<<"$line"
  for pair in "1:$R1" "2:$R2"; do
    m="${pair%%:*}"; f="${pair#*:}"
    [[ "$f" == - ]] && continue
    zcat "$f" | head -n $((4 * N)) | gzip -c > "$OUT/${S}_${m}.fq.gz" || true
    log "$S mate $m: $(( $(zcat "$OUT/${S}_${m}.fq.gz" | wc -l) / 4 )) reads"
  done
done
echo
echo "Add to config/datasets.tsv (tab-separated):"
printf 'TEST_%s\tTEST\tTEST_%s\t%s\t%s\tfirst %s reads of %s from %s\n' \
  "$SRC" "$SRC" "$(dataset_field "$SRC" layout)" "$(dataset_field "$SRC" strand)" "$N" "$*" "$SRC"
