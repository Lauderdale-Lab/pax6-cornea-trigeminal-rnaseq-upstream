#!/usr/bin/env bash
#
# tests/run_tests.sh -- tests that need no cluster, no data and no modules.
#
#   bash tests/run_tests.sh
#
# Stand-ins for sbatch, module, samtools and featureCounts (tests/stubs/) let
# the real code paths run: config parsing, sample discovery, the run and
# commit guards, job chaining, the strand gate before counting, the matrix
# tidy-up, and the count comparison. Run before every commit; CI runs it on
# every push.

set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PASS=0; FAIL=0
ok()   { PASS=$((PASS + 1)); printf '  ok    %s\n' "$1"; }
bad()  { FAIL=$((FAIL + 1)); printf '  FAIL  %s\n' "$1"; [[ -n "${2:-}" ]] && sed 's/^/          /' <<<"$2"; }
expect_ok()   { local d="$1"; shift; local o; if o=$("$@" 2>&1); then ok "$d"; else bad "$d" "$o"; fi; }
expect_fail() { local d="$1" pat="$2"; shift 2; local o; if o=$("$@" 2>&1); then bad "$d (should have failed)" "$o"
                elif grep -qE "$pat" <<<"$o"; then ok "$d"; else bad "$d (wrong message)" "$o"; fi; }

T=$(mktemp -d); [[ -n "${KEEP_TEST_DIR:-}" ]] && echo "test dir: $T" || trap 'rm -rf "$T"' EXIT
# Work on a committed copy of the repository so the git guards are real.
cp -r "$REPO" "$T/repo"; rm -rf "$T/repo/.git"
git -C "$T/repo" init -q; git -C "$T/repo" add -A
git -C "$T/repo" -c user.name=t -c user.email=t@t commit -qm init
P="$T/repo/bin/pax6"
export PAX6_ROOT="$T/root" PAX6_RUNS="$T/scratch/runs" PAX6_KEEP="$T/root/runs" PAX6_MAIL="" PATH="$REPO/tests/stubs:$PATH" SBATCH_LOG="$T/sbatch.log"
module() { return 0; }; export -f module

gz() { mkdir -p "$(dirname "$1")"; printf '@r\nACGT\n+\nIIII\n' | gzip -c > "$1"; }

echo "== syntax"
for f in "$REPO"/bin/pax6 "$REPO"/lib/*.sh "$REPO"/slurm/*.sbatch "$REPO"/tools/*.sh; do
  if bash -n "$f" 2>/dev/null; then ok "bash -n $(basename "$f")"; else bad "bash -n $f"; fi
done

echo "== config tables are read by column name"
cat > "$T/t.tsv" <<'EOF'
# comment
strand	dataset	layout
2	X	PE

0	Y	SE
EOF
(
  source "$REPO/lib/common.sh"
  [[ "$(tsv_col "$T/t.tsv" layout)" == 3 ]] && [[ "$(tsv_rows "$T/t.tsv" | wc -l)" == 2 ]] &&
  [[ "$(tsv_get "$T/t.tsv" 2 dataset)" == X ]]
) && ok "columns found by name, comments and blanks skipped" || bad "tsv helpers"

echo "== sample discovery"
R="$PAX6_ROOT/data/Lauderdale_2026/raw"
gz "$R/A_CW3/A_CW3_1.fq.gz"; gz "$R/A_CW3/A_CW3_2.fq.gz"
gz "$R/A_CU4/A_CU4_1.fq.gz"; gz "$R/A_CU4/A_CU4_2.fq.gz"
out=$("$P" samples Lauderdale_2026 2>&1)
if grep -q "A_CU4" <<<"$out" && grep -q "A_CW3/A_CW3_2.fq.gz" <<<"$out"; then ok "nested per-sample folders, mates paired"; else bad "nested discovery" "$out"; fi
gz "$PAX6_ROOT/data/Duncan_GSE183742/raw/SRR15826449_1.fastq.gz"; gz "$PAX6_ROOT/data/Duncan_GSE183742/raw/SRR15826449_2.fastq.gz"
out=$("$P" samples Duncan_GSE183742 2>&1)
grep -q "^SRR15826449 " <<<"$out" && ok "SRA-style .fastq.gz names" || bad "SRA names" "$out"

gz "$PAX6_ROOT/data/Lauderdale_2025/raw/B_CW1_1.fq.gz"
expect_fail "missing mate is refused" "No _2 mate" "$P" samples Lauderdale_2025
gz "$PAX6_ROOT/data/Lauderdale_2025/raw/B_CW1_2.fq.gz"
gz "$PAX6_ROOT/data/Lauderdale_2025/raw/copy/B_CW1_1.fq.gz"; gz "$PAX6_ROOT/data/Lauderdale_2025/raw/copy/B_CW1_2.fq.gz"
expect_fail "duplicate sample name is refused" "more than once" "$P" samples Lauderdale_2025
rm -rf "$PAX6_ROOT/data/Lauderdale_2025/raw/copy"
gz "$PAX6_ROOT/data/Lauderdale_2025/raw/stray.fq.gz"
expect_fail "mixed paired/single-end is refused" "mixes paired-end" "$P" samples Lauderdale_2025
rm -f "$PAX6_ROOT/data/Lauderdale_2025/raw/stray.fq.gz"

mkdir -p "$PAX6_ROOT/data/Lauderdale_2025/metadata"
printf 'sample_id\tr1\tr2\nB_CW1\tB_CW1_1.fq.gz\tB_CW1_2.fq.gz\n' > "$PAX6_ROOT/data/Lauderdale_2025/metadata/samples.tsv"
out=$("$P" samples Lauderdale_2025 2>&1)
grep -q "raw/B_CW1_1.fq.gz" <<<"$out" && ok "metadata/samples.tsv overrides discovery" || bad "manual samples" "$out"

echo "== strand call (fractions measured on the real libraries)"
(
  source "$REPO/lib/common.sh"
  [[ "$(strand_call 0.653 0.355 0.345)" == 0 ]] || { echo "Lauderdale A_CW1"; exit 1; }
  [[ "$(strand_call 0.674 0.350 0.352)" == 0 ]] || { echo "Lauderdale A_CW3"; exit 1; }
  [[ "$(strand_call 0.360 0.020 0.340)" == 2 ]] || { echo "Duncan"; exit 1; }
  [[ "$(strand_call 0.700 0.650 0.040)" == 1 ]] || { echo "forward"; exit 1; }
  [[ "$(strand_call 0.700 0.100 0.100)" == AMBIGUOUS ]] || { echo "ambiguous"; exit 1; }
) && ok "unstranded, reverse, forward and ambiguous patterns" || bad "strand_call"

echo "== runs and the commit guard"
expect_ok   "new-run creates a run"            "$P" new-run R1
[[ -f "$PAX6_RUNS/R1/RUN_INFO.tsv" && -f "$PAX6_RUNS/R1/modules.tsv" ]] && ok "run freezes RUN_INFO, modules, parameters" || bad "run files"
expect_fail "run names cannot be reused"        "already exists" "$P" new-run R1
expect_fail "unknown assembly is refused"       "Unknown assembly" "$P" new-run R2 hg38

echo "== submit chains the jobs"
: > "$SBATCH_LOG"
expect_ok "submit two datasets" "$P" submit R1 Lauderdale_2026 Duncan_GSE183742
n=$(wc -l < "$SBATCH_LOG")
[[ "$n" == 8 ]] && ok "4 jobs per dataset" || bad "expected 8 sbatch calls, got $n" "$(cat "$SBATCH_LOG")"
grep -q -- '--array=1-2%8' "$SBATCH_LOG" && ok "array sized to the libraries" || bad "array size" "$(cat "$SBATCH_LOG")"
grep 'align.sbatch' "$SBATCH_LOG" | grep -q -- '--dependency=aftercorr:' && ok "align waits per library on trim (aftercorr)" || bad "align dependency"
grep 'strand.sbatch' "$SBATCH_LOG" | grep -q -- '--dependency=afterok:' && ok "strand check waits on alignment" || bad "strand dependency"
grep -q 'PAX6_RUN=R1' "$SBATCH_LOG" && grep -q "PAX6_HOME=$T/repo" "$SBATCH_LOG" && ok "jobs receive run and repository" || bad "export"

gz "$R/A_CO5/A_CO5_1.fq.gz"; gz "$R/A_CO5/A_CO5_2.fq.gz"
expect_fail "new libraries cannot join a started run" "differ from those already" "$P" submit R1 Lauderdale_2026
rm -rf "$R/A_CO5"

echo "# a code change" >> "$T/repo/lib/common.sh"
expect_fail "uncommitted code change blocks submit" "Uncommitted changes" "$P" submit R1 Duncan_GSE183742
git -C "$T/repo" -c user.name=t -c user.email=t@t commit -qam change
expect_fail "committed code change blocks the old run" "changed since run R1" "$P" submit R1 Duncan_GSE183742
git -C "$T/repo" -c user.name=t -c user.email=t@t revert --no-edit HEAD >/dev/null
printf 'Extra\tExtra\tE\tPE\t0\tdoc\n' >> "$T/repo/config/datasets.tsv"
git -C "$T/repo" -c user.name=t -c user.email=t@t commit -qam data
: > "$SBATCH_LOG"
expect_ok "datasets.tsv edits are allowed mid-run" "$P" submit R1 Duncan_GSE183742

echo "== trim, align, strand check and MultiQC jobs, end to end on stand-in tools"
export EBROOTTRIMMOMATIC="$T/trimmomatic" TRIM_LOG="$T/trim.log"
mkdir -p "$EBROOTTRIMMOMATIC/adapters"; echo jar > "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"
G="$PAX6_ROOT/reference/GRCm38"; mkdir -p "$G/genome_snp_tran"
printf '>chr1\nACGT\n' > "$G/GRCm38.primary_assembly.genome.fa"
printf 'chr1\tHAVANA\texon\t1\t4\t.\t+\t.\tgene_id "g"; transcript_id "t";\n' > "$G/gencode.vM25.annotation.gtf"
echo idx > "$G/genome_snp_tran/genome_snp_tran.1.ht2"
printf '>a\nAGATCGGAAGAGC\n' > "$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
J() { env SLURM_ARRAY_TASK_ID="$1" PAX6_HOME="$T/repo" PAX6_RUN=R1 PAX6_DATASET=Lauderdale_2026 bash "$T/repo/slurm/$2"; }
D26="$PAX6_RUNS/R1/Lauderdale_2026"
expect_ok "trim task 1" J 1 qc_trim.sbatch
FASTQC_STUB_ADAPTER=1 expect_ok "trim task 2 (adapter present)" J 2 qc_trim.sbatch
[[ -s "$D26/trimmed/A_CU4_1.fq.gz" && -s "$D26/trimmed/A_CU4_2.fq.gz" && -f "$D26/trimmed/A_CU4.done" ]] && ok "trimmed pair and completion marker written" || bad "trim outputs" "$(ls -R "$D26")"
[[ -z "$(find "$D26/trimmed" -name '*unpaired*' -o -name '*.tmp.gz')" ]] && ok "no unpaired or temporary files left" || bad "leftovers" "$(ls "$D26/trimmed")"
[[ "$(ls "$D26/qc/trimmed"/*_fastqc.zip | wc -l)" == 4 ]] && ok "post-trim FastQC on both mates" || bad "post-trim FastQC" "$(ls "$D26/qc/trimmed")"
grep -q $'TruSeq3_default' "$D26/qc/adapters/A_CU4_adapter_source.tsv" && grep -q $'detected' "$D26/qc/adapters/A_CW3_adapter_source.tsv" && ok "adapter source recorded (default vs detected)" || bad "adapter source"
grep -q ':2:30:10:2:True$' "$TRIM_LOG" && ok "published ILLUMINACLIP settings used" || bad "clip settings" "$(cat "$TRIM_LOG")"
out=$(J 1 qc_trim.sbatch 2>&1); grep -q "already trimmed and verified" <<<"$out" && ok "resubmitted trim task skips verified output" || bad "trim resume" "$out"
# The two tasks start together, as array tasks do on the cluster. They once
# raced on one temporary splice-site file and the second failed.
J 1 align.sbatch > "$T/a1.log" 2>&1 & p1=$!
J 2 align.sbatch > "$T/a2.log" 2>&1 & p2=$!
if wait "$p1"; then ok "align task 1 (concurrent)"; else bad "align task 1" "$(cat "$T/a1.log")"; fi
if wait "$p2"; then ok "align task 2 (concurrent)"; else bad "align task 2" "$(cat "$T/a2.log")"; fi
[[ -z "$(find "$PAX6_RUNS/R1/reference" -name '*.tmp*')" ]] && ok "no temporary splice-site files left" || bad "splice tmp leftovers"
[[ -s "$D26/align/A_CU4.bam" && -s "$D26/align/A_CU4.bam.bai" && ! -e "$D26/align/A_CU4.bam.tmp" ]] && ok "BAM + index written under final name" || bad "align outputs" "$(ls "$D26/align")"
[[ -s "$PAX6_RUNS/R1/reference/splice_sites.txt" ]] && ok "splice sites extracted once per run" || bad "splice sites"
out=$(J 1 align.sbatch 2>&1); grep -q "verified; nothing to do" <<<"$out" && ok "resubmitted align task skips verified BAM" || bad "align resume" "$out"
out=$(env PAX6_HOME="$T/repo" PAX6_RUN=R1 PAX6_DATASET=Lauderdale_2026 bash "$T/repo/slurm/strand.sbatch" 2>&1)
grep -q "MEASURED   layout=PE strand=0" <<<"$out" && grep -q "^AGREE" <<<"$out" && ok "strand check measures unstranded PE and agrees with config" || bad "strand job" "$out"
expect_ok "MultiQC job" env PAX6_HOME="$T/repo" PAX6_RUN=R1 PAX6_DATASET=Lauderdale_2026 bash "$T/repo/slurm/multiqc.sbatch"
out=$("$P" status R1 2>&1)
grep -qE "^Lauderdale_2026 +2 +2 +2 +PE 0 +yes" <<<"$out" && ok "status reports libraries, trimmed, aligned, strand, MultiQC" || bad "status" "$out"

echo "== counting gates and outputs"
RUNR="$PAX6_RUNS/R1"
mkdir -p "$RUNR/Lauderdale_2026/align"
for s in A_CU4 A_CW3; do echo x > "$RUNR/Lauderdale_2026/align/$s.bam"; echo x > "$RUNR/Lauderdale_2026/align/$s.bam.bai"; done
mkdir -p "$RUNR/Lauderdale_2026/strand"; printf 'layout\tstrand\nPE\t2\n' > "$RUNR/Lauderdale_2026/strand/RESULT.tsv"
COUNT=(env PAX6_HOME="$T/repo" PAX6_RUN=R1 PAX6_STUDY=Lauderdale bash "$T/repo/slurm/count.sbatch")
expect_fail "no strand result for a study dataset" "no strand result" "${COUNT[@]}"
# Lauderdale also has 2025; give it a run directory with one library
mkdir -p "$RUNR/Lauderdale_2025/align" "$RUNR/Lauderdale_2025/strand"
printf 'B_CW1\t-\t-\n' > "$RUNR/Lauderdale_2025/samples.tsv"
echo x > "$RUNR/Lauderdale_2025/align/B_CW1.bam"; echo x > "$RUNR/Lauderdale_2025/align/B_CW1.bam.bai"
printf 'layout\tstrand\nPE\t0\n' > "$RUNR/Lauderdale_2025/strand/RESULT.tsv"
expect_fail "measured strand disagreeing with config stops counting" "measured layout=PE strand=2" "${COUNT[@]}"
printf 'layout\tstrand\tassigned_s0\tassigned_s1\tassigned_s2\tlibraries_tested\nPE\t0\t0.674\t0.350\t0.352\t2\n' > "$RUNR/Lauderdale_2026/strand/RESULT.tsv"
expect_ok "counting runs when every dataset agrees" "${COUNT[@]}"
C="$RUNR/counts/Lauderdale"
[[ "$(head -1 "$C/gene_counts.tsv")" == $'gene_id\tB_CW1\tA_CU4\tA_CW3' ]] && ok "tidy matrix: sample IDs as columns, in config order" || bad "matrix header" "$(head -2 "$C/gene_counts.tsv")"
grep -q $'^strand\t0$' "$C/PROVENANCE.tsv" && grep -q '^matrix_md5' "$C/PROVENANCE.tsv" && ok "provenance records strand and checksum" || bad "provenance" "$(cat "$C/PROVENANCE.tsv")"
grep -q -- '--countReadPairs' "$C/PROVENANCE.tsv" && grep -q -- '-s 0' "$C/PROVENANCE.tsv" && ok "fragments counted at the measured strand" || bad "command" "$(cat "$C/PROVENANCE.tsv")"
PV="$C/provenance"
for f in software_versions.tsv reference.tsv inputs.tsv library_qc.tsv strand.tsv METHODS.md parameters.env modules.tsv RUN_INFO.tsv; do
  [[ -s "$PV/$f" ]] || bad "provenance/$f missing"
done
grep -qP '^HISAT2\tHISAT2/2.2.1-gompi-2023a\thisat2-align-s version 2.2.1\talign' "$PV/software_versions.tsv" &&
  grep -qP '^Trimmomatic\t.*\tTrimmomatic 0.39\t' "$PV/software_versions.tsv" &&
  grep -qP '^FastQC\t.*\tFastQC v0.12.1\t' "$PV/software_versions.tsv" &&
  grep -qP '^Subread\t.*\tfeatureCounts v2.0.6\tcount,strand|^Subread\t.*\tfeatureCounts v2.0.6\tstrand,count' "$PV/software_versions.tsv" &&
  ok "versions as the tools report them, with the steps that used them" || bad "software_versions" "$(cat "$PV/software_versions.tsv")"
grep -qP '^annotation\t.*\t[0-9a-f]{32}\t[0-9]+ genes' "$PV/reference.tsv" && grep -qP '^index_file\t.*\.1\.ht2\t[0-9a-f]{32}' "$PV/reference.tsv" &&
  ok "reference files checksummed" || bad "reference.tsv" "$(cat "$PV/reference.tsv")"
grep -qP '^A_CW3\tLauderdale_2026\t2026\t.*A_CW3_1.fq.gz\t[0-9a-f]{32}\t.*A_CW3_2.fq.gz\t[0-9a-f]{32}$' "$PV/inputs.tsv" &&
  ok "raw FASTQ checksums per library" || bad "inputs.tsv" "$(cat "$PV/inputs.tsv")"
grep -qP '^A_CW3\tLauderdale_2026\t2026\t1\t100.00\tdetected\t97.00\t90.00\t900\t' "$PV/library_qc.tsv" &&
  grep -qP '^B_CW1\tLauderdale_2025\t2025\tNA\t' "$PV/library_qc.tsv" &&
  ok "per-library QC table (NA where a log is missing, not a crash)" || bad "library_qc.tsv" "$(cat "$PV/library_qc.tsv")"
grep -q "aligned to the GRCm38 primary assembly" "$PV/METHODS.md" && grep -q "with HISAT2 2.2.1 (--dta)" "$PV/METHODS.md" && grep -q "featureCounts (Subread 2.0.6)" "$PV/METHODS.md" &&
  grep -q "transcript-aware HISAT2 index" "$PV/METHODS.md" && grep -q "were unstranded (featureCounts -s 0)" "$PV/METHODS.md" &&
  grep -q "from 3 libraries" "$PV/METHODS.md" && ! grep -q "not recorded" "$PV/METHODS.md" &&
  ok "Methods draft filled in from the record" || bad "METHODS.md" "$(cat "$PV/METHODS.md")"
K="$PAX6_KEEP/R1"
[[ -s "$K/counts/Lauderdale/gene_counts.tsv" && -s "$K/counts/Lauderdale/provenance/METHODS.md" && -s "$K/RUN_INFO.tsv" ]] &&
  ok "matrix, provenance and run records kept on /work automatically" || bad "keep after count" "$(find "$PAX6_KEEP" -type f 2>/dev/null | head)"
[[ -z "$(find "$K" -name '*.bam' -o -path '*/trimmed/*' -type f)" ]] && ok "BAMs and trimmed reads are not copied by default" || bad "keep copied BAMs or trimmed reads"
cmp -s "$C/gene_counts.tsv" "$K/counts/Lauderdale/gene_counts.tsv" && ok "kept matrix identical to the one on scratch" || bad "kept matrix differs"
expect_ok "keep_run.sh --with-bams" "$T/repo/tools/keep_run.sh" R1 --with-bams
[[ -s "$K/Lauderdale_2026/align/A_CW3.bam" ]] && ok "keep_run.sh --with-bams copies BAMs" || bad "BAMs not kept"
expect_fail "an existing matrix is never overwritten" "never overwritten" "${COUNT[@]}"
rm -rf "$C"; expect_fail "low assignment ratio fails the job" "below 70%" env FC_STUB_LOW=1 "${COUNT[@]}"

echo "== compare_counts.py"
printf '# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t/x/A_CW3_Novogene_X202SC26055473_hisat2.sorted.bam\t/y/B_CW1_hisat2.sorted.bam\nENSMUSG01.3\t1\t1\t2\t+\t2\t10\t5\nENSMUSG02.1\t1\t1\t2\t+\t2\t0\t7\n' > "$T/old.txt"
printf 'gene_id\tA_CW3\tB_CW1\nENSMUSG01.4\t10\t6\nENSMUSG02.1\t0\t7\n' > "$T/new.tsv"
out=$(python3 "$REPO/tools/compare_counts.py" "$T/new.tsv" "$T/old.txt" 2>&1)
if grep -q "2 in both" <<<"$out" && grep -q "1 of 2 libraries identical" <<<"$out"; then ok "old and new naming line up; differences found"; else bad "compare_counts" "$out"; fi

echo
echo "$PASS passed, $FAIL failed"
(( FAIL == 0 ))
