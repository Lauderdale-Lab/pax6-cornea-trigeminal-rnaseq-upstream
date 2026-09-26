#!/usr/bin/env bash
# Project audit. Read-only: reports, changes nothing.
R=/work/jdllab/PAX6_RNAseq
ok(){ printf '  OK    %s\n' "$*"; }
bad(){ printf '  MISS  %s\n' "$*"; }
chk(){ [ -e "$2" ] && ok "$1" || bad "$1  ($2)"; }
n(){ printf '%s' "$(find -L "$1" -maxdepth "${3:-2}" -type f -name "$2" 2>/dev/null | wc -l)"; }

echo "=== 1. top level ==="
for d in pipeline docs data counts reference logs; do chk "$d/" "$R/$d"; done

echo; echo "=== 2. pipeline: 10 files ==="
for f in pipeline_common.sh studies.tsv assemblies.tsv module_versions.txt \
         build_hisat2_index.sh qc_trim_batch.sh \
         align_batch.sh strand_check.sh count_study.sh sam_to_bam.sh; do
  chk "$f" "$R/pipeline/$f"
done
echo "  stray logs in pipeline/: $(ls $R/pipeline/*.out $R/pipeline/*.err 2>/dev/null | wc -l) (want 0)"
echo "  syntax:"; for f in $R/pipeline/*.sh; do bash -n "$f" 2>/dev/null || echo "    FAIL $(basename $f)"; done

echo; echo "=== 3. manifest resolves ==="
awk -F'\t' '$0!~/^#/ && NF>=6 {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $R/pipeline/studies.tsv |
while IFS=$'\t' read -r st ba lay str tr map; do
  printf '  %-11s %-10s layout=%-4s strand=%-4s ' "$st" "$ba" "$lay" "$str"
  [ -d "$tr" ] && printf 'trimmed=OK ' || printf 'trimmed=MISS '
  [ -d "$map/GRCm38" ] && printf 'bams=OK\n' || printf 'bams=MISS\n'
done

echo; echo "=== 4. datasets ==="
for D in $R/data/*/; do
  b=$(basename "$D")
  printf '  %-18s prov=%-4s raw=%-3s trim=%-4s bam=%-3s meta=%-3s ro=%s\n' "$b" \
    "$( [ -f "$D/PROVENANCE.md" ] && wc -l < "$D/PROVENANCE.md" || echo NO )" \
    "$(n "$D/raw" '*.f*q.gz' 3)" "$(n "$D/trimmed" '*.f*q.gz' 2)" \
    "$(n "$D/align/GRCm38" '*.sorted.bam' 1)" "$(ls "$D/metadata" 2>/dev/null | wc -l)" \
    "$( [ -w "$D/raw" ] && echo WRITABLE || echo readonly )"
done

echo; echo "=== 5. counts ==="
for C in $R/counts/*/; do
  printf '  %-30s ' "$(basename $C)"
  m="$C/gene_counts_featureCounts.txt"
  [ -f "$m" ] && printf 'matrix=%s cols=%s ' "$(du -h "$m"|cut -f1)" "$(sed -n '2p' "$m"|awk -F'\t' '{print NF-6}')" || printf 'NO MATRIX '
  [ -f "$C/sample_map.tsv" ] && echo "map=OK" || echo "map=MISS"
done

echo; echo "=== 6. reference ==="
for A in GRCm38 GRCm39; do
  printf '  %-8s fasta=%s gtf=%s ht2=%s\n' "$A" \
    "$(ls $R/reference/$A/*.fa 2>/dev/null | wc -l)" \
    "$(ls $R/reference/$A/*.gtf 2>/dev/null | wc -l)" \
    "$(ls $R/reference/$A/*/*.ht2 2>/dev/null | wc -l)"
done

echo; echo "=== 7. docs ==="
for f in PROJECT_LAYOUT.md ADDING_A_BATCH.md; do chk "$f" "$R/docs/$f"; done

echo; echo "=== 8. leftovers under /work/jdllab ==="
du -sh /work/jdllab/* 2>/dev/null | sort -h | tail -8
