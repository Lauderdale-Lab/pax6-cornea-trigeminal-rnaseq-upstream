#!/usr/bin/env bash
#
# tools/verify_md5.sh <dataset> -- check raw FASTQ against the vendor's MD5
# manifests (any MD5.txt, md5.txt or *.md5 under data/<dataset>/raw/).
# Run as a job for large datasets:  sbatch --mem=2G --time=8:00:00 tools/verify_md5.sh <dataset>
#
#SBATCH --job-name=pax6_md5
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00

set -uo pipefail
PAX6_HOME="${PAX6_HOME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source=../config/site.env
source "$PAX6_HOME/config/site.env"
DS="${1:?usage: verify_md5.sh <dataset>}"
RAW="$PAX6_ROOT/data/$DS/raw"
[[ -d "$RAW" ]] || { echo "No raw directory: $RAW" >&2; exit 1; }

ok=0; bad=0; n=0
while IFS= read -r m; do
  n=$((n + 1))
  if out=$(cd "$(dirname "$m")" && md5sum -c "$(basename "$m")" 2>&1); then
    ok=$((ok + 1)); echo "OK    ${m#"$RAW"/}"
  else
    bad=$((bad + 1)); echo "FAIL  ${m#"$RAW"/}"; sed 's/^/        /' <<<"$out"
  fi
done < <(find -L "$RAW" \( -iname 'md5.txt' -o -name '*.md5' \) -type f | sort)

(( n > 0 )) || { echo "No MD5 manifests found under $RAW (GEO/SRA downloads usually have none)."; exit 2; }
echo; echo "$DS: $ok manifest(s) OK, $bad failed"
(( bad == 0 ))
