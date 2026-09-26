# Adding a batch or a study

Study: *Nerve remodeling in a Pax6 model of keratopathy*
Last revised: 25 September 2026

The only file you edit is `studies.tsv`. Add one row, then run four jobs.

```bash
cd /work/jdllab/PAX6_RNAseq/pipeline

# Row format, tab separated:
#   study <TAB> batch <TAB> auto <TAB> auto <TAB> trimmed_dir <TAB> mapping_root
# layout and strand start as `auto`; strand_check.sh fills them in.

sbatch qc_trim_batch.sh  Lauderdale 2027 /path/to/raw   # QC + trim + MultiQC
sbatch align_batch.sh    Lauderdale 2027                # BAMs; array-capable
sbatch strand_check.sh   Lauderdale 2027                # reports layout AND strand
# record layout and strand in studies.tsv from the RESULT line
sbatch count_study.sh    Lauderdale                     # whole study, one matrix
```

Every script takes `<study> <batch>` — both, always. `count_study.sh` takes the
study alone, because counting spans all of that study's batches.

**Batch or study?** A new batch of your own libraries — same alleles, same
design — is a new *batch* of the `Lauderdale` study, and counting re-runs over
all of its batches. A public dataset, or anything carrying a different allele,
is a new *study*: quantified identically so results are comparable, counted and
modelled separately. `count_study.sh` never mixes them.

**Which assembly?** Everything defaults to the first row of `assemblies.tsv`
(GRCm38 + GENCODE vM25, the build the published results use). For another:

```bash
sbatch --export=ALL,PAX6_ASSEMBLY=GRCm39 align_batch.sh Lauderdale 2026
sbatch --export=ALL,PAX6_ASSEMBLY=GRCm39 count_study.sh Lauderdale
```

BAMs go to `<mapping_root>/<assembly>/` and the count directory carries the
assembly in its name, so the same libraries can exist against both builds at
once and a matrix can never contain a mixture. A new assembly needs a one-time
`sbatch build_hisat2_index.sh <assembly>` after its genome and matching
annotation are in `reference/<assembly>/`.

**SAM instead of FASTQ?** `sam_to_bam.sh <study> <batch> <sam_dir>` converts,
but refuses to run until you have inspected the `@PG` line it prints and
confirmed the alignment used this project's reference. Converting inherits an
unknown alignment; re-aligning from FASTQ is what you want for anything you
intend to publish.

Afterwards, follow the "NEXT" block the counting job prints: repoint
`FILE_COUNTS`, archive the superseded matrix, add the new samples to the
metadata CSV **with their batch label**, and re-run the analysis into a new
output subdirectory.

---

## Running it on Sapelo2

### One-time setup

Ten files live in `/work/jdllab/PAX6_RNAseq/pipeline` — the seven scripts in
this repository's `pipeline/` folder plus `studies.tsv`, `assemblies.tsv` and
`module_versions.txt`. Keep them on `/work`, not `/home` (quota) and not
`/scratch` (purged). The simplest way to keep the cluster copy identical to
this repository is to `git clone` it there and submit from `pipeline/`.

If they came from macOS, strip carriage returns once — a stray `\r` makes bash
fail with unreadable errors:

```bash
cd /work/jdllab/PAX6_RNAseq/pipeline
sed -i 's/\r$//' *.sh *.tsv module_versions.txt
chmod +x *.sh
```

**Always `cd` to the pipeline directory before submitting.** The scripts locate
each other through `SLURM_SUBMIT_DIR`, and job logs land in the submit
directory.

### Check before you queue

All of this runs on the login node in seconds. Never run alignment or counting
there.

```bash
cd /work/jdllab/PAX6_RNAseq/pipeline

# 1. Syntax.
for f in *.sh; do bash -n "$f" || echo "SYNTAX: $f"; done

# 2. Every manifest path.
awk -F'\t' '$0 !~ /^#/ && NF>=6 { print $1, $2, $5, $6 }' studies.tsv |
  while read -r study batch trimmed mapping; do
    [[ -d "$trimmed" ]] && echo "OK   $study/$batch trimmed" \
                        || echo "MISS $study/$batch trimmed  $trimmed"
    [[ -d "$mapping/GRCm38" ]] && echo "OK   $study/$batch BAMs (GRCm38)" \
                        || echo "note $study/$batch BAM dir will be created: $mapping/GRCm38"
  done

# 3. Pinned modules still on the cluster.
while read -r tool ver; do
  [[ "$tool" == \#* || -z "$tool" ]] && continue
  ml "$tool/$ver" 2>/dev/null && echo "OK   $tool/$ver" || echo "MISS $tool/$ver"
done < module_versions.txt
ml purge

# 4. Reference files match assemblies.tsv.
ls /work/jdllab/PAX6_RNAseq/reference/GRCm38/
ls /work/jdllab/PAX6_RNAseq/reference/GRCm38/genome_snp_tran/*.1.ht2
```

### Submitting

```bash
cd /work/jdllab/PAX6_RNAseq/pipeline

sbatch align_batch.sh Lauderdale 2027                  # whole batch, one job
sbatch --array=1-21%6 align_batch.sh Lauderdale 2027   # one sample per task
```

Use array mode for anything beyond a handful of libraries: a single failure
costs one sample rather than the batch, and samples run in parallel. Set the
range to the library count; `%6` caps concurrency.

Chaining alignment to the strand check is safe — the strand check only needs
BAMs:

```bash
JID=$(sbatch --parsable align_batch.sh Lauderdale 2027)
sbatch --dependency=afterok:$JID strand_check.sh Lauderdale 2027
```

Do **not** chain the count job. It needs a strand value a human has read and
recorded in `studies.tsv`.

### Monitoring

```bash
squeue --me
sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
seff <jobid>                   # after it finishes: memory and CPU actually used
tail -f align_<jobid>_*.out    # live progress
```

Check `seff` on the first run of a new batch. If `MaxRSS` sits far below the
48G requested, lower it and jobs start sooner; if a job was killed for memory,
`sacct` shows `OUT_OF_MEMORY` rather than a useful error.

### Resource notes

Sized for these libraries (50–70M pairs, 2 × 150 bp):

| Job | cpus | mem | walltime | actual |
|---|---|---|---|---|
| `qc_trim_batch.sh` | 8 | 32G | 48h | ~1 h per library |
| `align_batch.sh` | 12 | 48G | 48h | ~2–3 h per library |
| `strand_check.sh` | 8 | 16G | 2h | ~5 min |
| `count_study.sh` | 8 | 32G | 8h | ~40 min for 48 libraries |
| `build_hisat2_index.sh` | 16 | 64G | 72h | hours (plain index) |

The 2026 batch of 21 libraries took about 2.5 days as a single serial job. As
an array it would have taken an afternoon.

If a job hits the walltime nothing is lost — resubmitting skips every verified
BAM and resumes.

---

## Why the pipeline is shaped this way

Each guard exists because the corresponding thing went wrong.

**Module versions are pinned** (`module_versions.txt`). The old scripts
hard-coded some modules, loaded "the latest available version" of others, and
recorded neither in their logs. At Subread 2.0.2 the meaning of `-p` changed
from "count fragments" to "the input is paired-end". The old counting command
was written for the earlier meaning and ran under Subread 2.0.6 in both
batches, so it produced a matrix of reads rather than fragments, in which `-B`
and `-C` never took effect. Nothing failed, and it went unnoticed for months.
A pinned version that is reported in every log, and fails loudly when
unavailable, makes that kind of change visible.

**Strandedness and layout are measured, never assumed.** Both existing batches
were counted at `-s 2` on the assumption of a directional dUTP library. They are
unstranded. Roughly half of every library's exon-overlapping fragments were
discarded. The error was invisible downstream precisely because it applied to
every sample equally — it distorted no comparison, it only destroyed power.
`strand_check.sh` costs minutes; skipping it cost this study a full re-count.

**Counting spans batches by construction.** It used to live inside the
per-batch alignment script, so the two batches were counted months apart in
separate jobs — a processing difference that could sit exactly on top of
the batch covariate that appears in every model. Adding a batch now means
re-counting the whole study, which is the point.

**Resume checks the product, not the filename.** The old scripts skipped a step
whenever an output path existed, so a BAM truncated by a killed job counted as
finished. `bam_ok()` requires `samtools quickcheck` to pass and an index to
exist.

**One naming convention, and batch lives in the manifest.** The 2026 BAMs
embedded a project tag in their filenames, which then had to be unpicked with
string surgery in R. New BAMs are `<sample_id>_hisat2.sorted.bam`;
`sample_map.tsv` carries study, batch and assembly. Counting aborts if a sample
ID appears in two batches of one study.

**The intermediate SAM is never written.** `hisat2` pipes straight into
`samtools sort`. At this depth each SAM ran to roughly 200 GB, and the old
pipeline kept all of them.

**Transcript assembly was removed.** StringTie and prepDE produced output no
analysis ever consumed, and both were mis-parameterised — `--rf` on unstranded
libraries, prepDE at its default 75 bp against 150 bp reads. Wrong files
sitting in an output directory are a liability, because a later script may find
them. (Note: "assembly" elsewhere in this pipeline means the *genome* build.)

---

## When something stops

Every failure below is deliberate. None should be worked around by loosening
the check.

| Message | What it means |
|---|---|
| `No pinned version for 'X'` | A tool has no entry in `module_versions.txt`. |
| `Module not available: X/Y` | The cluster retired that version. Decide whether to move the pin — and if you do, re-count **every** batch of the study. |
| `Unknown assembly 'X'` | Not a row in `assemblies.tsv`; the message lists what is. |
| `Study 'X' is not in studies.tsv` | Typo, or the row was never added. |
| `<study>/<batch> has strand=auto` | Run `strand_check.sh <study> <batch>` and record the answer. |
| `<study>/<batch> has layout=auto` | Same job reports it; record PE or SE. |
| `Batches of 'X' disagree on strandedness` | One invocation applies one `-s`. Do not pick one and proceed; count the groups separately and handle the difference explicitly. |
| `Batches of 'X' mix PE and SE` | Fragment counts and read counts are not on the same scale. Split the study. |
| `no BAMs for assembly X at <dir>` | That batch has not been aligned to that assembly yet. |
| `BAM incomplete, unindexed or truncated` | Re-run `align_batch.sh`; it redoes only what is missing. |
| `Sample IDs occur in more than one batch` | Rename before counting. Two columns mapping to one metadata row is silent corruption. |
| `RESULT: AMBIGUOUS` from strand_check | The three settings match no expected pattern. Inspect the `.summary` files. Do not guess. |
| `[STOP] Provenance not confirmed` | `sam_to_bam.sh` wants you to read the `@PG` line first. |
| Assignment ratio `<-- LOW` | Near 50% across the board means the strand setting is wrong. A single low library is more likely degraded or low-complexity RNA. |

## Expected values, for reference

From the 48-library Lauderdale run of 13 August 2026, GRCm38, unstranded (`-s 0`),
paired-end fragments:

- Assigned 65–78% (median ~70%)
- Assigned ÷ (Assigned + NoFeatures) 81–97% (median ~94%)
- `Unassigned_MappingQuality` ~10–12% — `-Q 10` removes most of what `-M`
  admits, since HISAT2 gives multi-mapping alignments low MAPQ, so this
  effectively counts near-uniquely-aligned fragments
- `Unassigned_Ambiguity` and `Unassigned_MultiMapping` both zero
- Overall alignment rate 94–98%
- Library sizes 29–59M fragments, 14,185 genes with no counts anywhere
