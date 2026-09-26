#!/usr/bin/env bash
#SBATCH --job-name=verify_md5
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=verify_md5_%j.out
#SBATCH --error=verify_md5_%j.err

set -uo pipefail
D="${1:-/work/jdllab/PAX6_RNAseq/data/Lauderdale_2026/raw}"
cd "$D" || exit 1
ok=0; bad=0; missing=0
for s in */; do
  s="${s%/}"
  if [[ ! -f "$s/MD5.txt" ]]; then echo "NO MANIFEST: $s"; missing=$((missing+1)); continue; fi
  out=$( cd "$s" && md5sum -c MD5.txt 2>&1 )
  if grep -q "FAILED\|WARNING" <<< "$out"; then
    echo "FAIL: $s"; echo "$out" | sed 's/^/    /'; bad=$((bad+1))
  else
    echo "OK:   $s"; ok=$((ok+1))
  fi
done
echo
echo "=== $ok samples OK, $bad failed, $missing without a manifest ==="
