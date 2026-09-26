#!/usr/bin/env bash
#
# tools/keep_run.sh <run> [--with-bams] -- copy a run's results from scratch
# to /work ($PAX6_KEEP/<run>). Safe to repeat: only new or changed files are
# copied.
#
# Count jobs already do this for everything except BAMs. Run it yourself to
# keep the BAMs (needed for any later per-allele or per-exon work on these
# libraries), or to refresh the QC and logs before scratch is purged.
#
#   tools/keep_run.sh reprocess_2026-09 --with-bams

set -euo pipefail
PAX6_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../lib/common.sh
source "$PAX6_HOME/lib/common.sh"
RUN="${1:?usage: keep_run.sh <run> [--with-bams]}"
dst=$(keep_results "$RUN" "${2:-}")
log "kept: $dst ($(du -sh "$dst" 2>/dev/null | cut -f1))"
