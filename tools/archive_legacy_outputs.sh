#!/usr/bin/env bash
#
# tools/archive_legacy_outputs.sh -- move everything the August 2026 pipeline
# produced out of the way, so $PAX6_ROOT holds only inputs (data/, reference/)
# and new-pipeline outputs (runs/).
#
#   tools/archive_legacy_outputs.sh            DRY RUN: prints, changes nothing
#   tools/archive_legacy_outputs.sh --apply    moves
#
# Moves (renames within /work, so instant and needing no extra space):
#   data/<dataset>/{trimmed,qc,align}   -> archive/2026-08_pipeline/data/<dataset>/
#   counts/                             -> archive/2026-08_pipeline/counts/
#   pipeline/ docs/ logs/               -> archive/2026-08_pipeline/
# Never touches data/<dataset>/raw, data/<dataset>/metadata, PROVENANCE.md or
# reference/. Deletes nothing. Refuses to overwrite anything in the archive.

set -euo pipefail
PAX6_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../config/site.env
source "$PAX6_HOME/config/site.env"
APPLY=0; [[ "${1:-}" == --apply ]] && APPLY=1
DEST="$PAX6_ROOT/archive/2026-08_pipeline"
MANIFEST="$DEST/ARCHIVE_MANIFEST.tsv"

(( APPLY )) || echo "=== DRY RUN -- nothing will change. Re-run with --apply. ==="
echo "project root: $PAX6_ROOT"
echo "archive     : $DEST"
echo

# This repository must not be moved while it is running.
case "$PAX6_HOME/" in "$PAX6_ROOT/pipeline/"*)
  echo "[FATAL] This copy of the pipeline lives in $PAX6_ROOT/pipeline, which is being archived." >&2
  echo "        Clone the new pipeline elsewhere first (e.g. $PAX6_ROOT/pipeline_new) and run it from there." >&2
  exit 1 ;;
esac

n=0
move() {
  local src="$1" dst="$2"
  [[ -e "$src" || -L "$src" ]] || return 0
  if [[ -e "$dst" || -L "$dst" ]]; then echo "  SKIP  $src  (already in archive: $dst)"; return 0; fi
  echo "  move  ${src#"$PAX6_ROOT"/}  ->  ${dst#"$PAX6_ROOT"/}"
  n=$((n + 1))
  if (( APPLY )); then
    mkdir -p "$(dirname "$dst")"
    mv "$src" "$dst"
    printf '%s\t%s\t%s\n' "$(date -Iseconds)" "$src" "$dst" >> "$MANIFEST"
  fi
}

(( APPLY )) && { mkdir -p "$DEST"; [[ -f "$MANIFEST" ]] || printf 'moved_at\tfrom\tto\n' > "$MANIFEST"; }

for D in "$PAX6_ROOT"/data/*/; do
  D="${D%/}"; ds=$(basename "$D")
  for sub in trimmed qc align; do move "$D/$sub" "$DEST/data/$ds/$sub"; done
done
for top in counts pipeline docs logs; do move "$PAX6_ROOT/$top" "$DEST/$top"; done

echo
if (( n == 0 )); then echo "nothing to move"; exit 0; fi
if (( APPLY )); then
  echo "$n item(s) moved; record in $MANIFEST"
  echo "The published matrices are now under $DEST/counts/ -- tools/compare_counts.py reads them from there."
else
  echo "$n item(s) would be moved. Re-run with --apply."
fi
