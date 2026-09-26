# Runbook — Sapelo2

Everything is run from the pipeline checkout with `bin/pax6`.

**Where things go.** Runs are processed on scratch,
`/scratch/$USER/PAX6_RNAseq/runs/<run>/`, where trimmed reads, BAMs and job
logs (`logs/`) are written. Scratch is fast, but Sapelo2 purges it. Every count
job therefore copies its matrix, the provenance record, and the run's QC and
records to `/work/jdllab/PAX6_RNAseq/runs/<run>/` the moment it finishes. To
keep the BAMs as well, run `tools/keep_run.sh <run> --with-bams`. Both
locations are set in `config/site.env`.

## 1. One-time setup: archive the old outputs, install the pipeline

The August 2026 pipeline left trimmed reads, QC, BAMs and count matrices
inside `data/`, plus old scripts in `pipeline/`. Archive them so `data/` holds
only inputs, then install this repository as `pipeline/`.

```bash
cd /work/jdllab/PAX6_RNAseq
git clone https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream.git pipeline_new
cd pipeline_new

tools/archive_legacy_outputs.sh            # dry run: read the list
tools/archive_legacy_outputs.sh --apply    # moves; deletes nothing

cd .. && mv pipeline_new pipeline && cd pipeline
bin/pax6 check
```

The archive lands in `archive/2026-08_pipeline/`, with a manifest of every
move. The published matrices end up under `archive/2026-08_pipeline/counts/`.
Moves within `/work` are renames: instant, and needing no extra space.

`check` should show all three datasets with their library counts, the GRCm38
reference files, and every pinned module as `OK`.

## 2. Confirm the sample IDs

```bash
bin/pax6 samples Lauderdale_2025
bin/pax6 samples Lauderdale_2026
bin/pax6 samples Duncan_GSE183742
```

Sample IDs are taken from file names (`<id>_1.fq.gz`). They become the
matrix's column names, so they must match the sample metadata the analyses
use (`A_CW3`, not `A_CW3_Novogene_X202SC26055473`).

If a dataset's file names don't give the right IDs, write
`data/<dataset>/metadata/samples.tsv` and it takes precedence over the file
names. The same applies if a library was sequenced across several lanes; merge
the lanes first with `cat`.

```
sample_id	r1	r2
A_CW3	A_CW3/A_CW3_EKRN2600_L1_1.fq.gz	A_CW3/A_CW3_EKRN2600_L1_2.fq.gz
```

Paths are relative to `data/<dataset>/raw/`.

## 3. Test on a small dataset first

This takes a few minutes of cluster time and checks every step on real reads.

```bash
tools/make_test_dataset.sh Lauderdale_2026 200000 A_CW3 A_CU3
# paste the printed row into config/datasets.tsv, then commit it
git commit -am "Add TEST_Lauderdale_2026 dataset"

bin/pax6 new-run test_$(date +%Y%m%d)
bin/pax6 submit  test_$(date +%Y%m%d) TEST_Lauderdale_2026
bin/pax6 status  test_$(date +%Y%m%d)          # repeat until strand and multiqc appear
bin/pax6 count   test_$(date +%Y%m%d) TEST
```

## 4. Reprocess everything

```bash
RUN=reprocess_$(date +%Y-%m)
bin/pax6 new-run $RUN
bin/pax6 submit  $RUN Lauderdale_2025 Lauderdale_2026 Duncan_GSE183742
bin/pax6 status  $RUN
```

`count` can be submitted at any time after `submit`: it waits for any strand
check of the study that is still queued or running. It stops cleanly if a
strand check failed, or if a measured strand disagrees with
`config/datasets.tsv`. Each strand log
(`/scratch/$USER/PAX6_RNAseq/runs/$RUN/logs/strand_<dataset>_*.out`) prints
`AGREE` or `DISAGREE`. To count:

```bash
bin/pax6 count $RUN Lauderdale
bin/pax6 count $RUN Duncan
```

Then compare with the published matrices before anyone uses the new numbers:

```bash
# needs only python3 (the system one, or: ml Python)
K=/work/jdllab/PAX6_RNAseq/runs/$RUN          # the kept copy on /work
tools/compare_counts.py $K/counts/Lauderdale/gene_counts.tsv \
  ../archive/2026-08_pipeline/counts/Lauderdale_GRCm38_20260813/gene_counts_featureCounts.txt \
  --out $K/counts/Lauderdale/compare_with_published.tsv
```

Identical trimming, alignment and counting settings should give nearly
identical counts. Any difference must be explained before the new matrix
replaces the old.

After counting, keep the BAMs on /work if later work will need them, such as
per-exon or per-allele analysis of *Pax6*. This is about 3 GB per library:

```bash
tools/keep_run.sh $RUN --with-bams
```

Check the space first with `lfs quota -h -g jdllab /work`. Once the new BAMs
are kept and the comparison is settled, the archived August BAMs
(`archive/2026-08_pipeline/data/*/align/`, about 160 GB) are the obvious space
to reclaim.

## 5. Adding a new dataset

1. Put its FASTQ in `data/<Lab>_<year|accession>/raw/` and make that folder
   read-only (`chmod -R a-w`) once `tools/verify_md5.sh <dataset>` passes.
2. Write `data/<dataset>/PROVENANCE.md` from `docs/PROVENANCE_template.md`.
3. Add a row to `config/datasets.tsv` with `layout` and `strand` set to
   `auto`, and commit it.
4. Process it in a **new** run: `new-run`, then `submit`.
5. Record the measured layout and strand in `config/datasets.tsv` and commit.
6. Run `count` for its study.

A new batch of the same Sey-Neu design joins the `Lauderdale` study; the whole
study is then re-counted in that run. A different allele or an outside lab
gets a study of its own.

## Resources

| Job | cpus | mem | time | Expect |
|---|---|---|---|---|
| qc_trim | 8 | 16G | 12h | ~1 h per library |
| align | 12 | 32G | 12h | ~1–2 h per library |
| strand | 8 | 8G | 2h | minutes |
| multiqc | 2 | 8G | 2h | minutes |
| count | 8 | 16G | 8h | ~40 min for 48 libraries |

Arrays run at most `PAX6_ARRAY_MAX` (8) libraries at once. After the first
real run, `seff <jobid>` shows memory actually used. Adjust the `#SBATCH`
lines, and commit the change, if the requests are far off.

If a task fails or hits its time limit, fix the cause and run `submit` again
with the same run and dataset. Verified outputs are skipped. Cancel any jobs
left waiting on the failed one first: `squeue --me`, then `scancel <jobid>`.

## Why it stops

Every stop is deliberate. Don't loosen a check to get past it.

| Message | Meaning | Guards against |
|---|---|---|
| `Uncommitted changes in pipeline code or parameters` | Commit before creating or continuing a run | a run that no commit describes |
| `changed since run <run> was created` | Code or parameters differ from the run's commit | one run processed by two pipeline versions |
| `libraries found differ from those already in this run` | Files appeared or vanished in `raw/` | a run whose sample set changed partway |
| `No pinned module` / `Module not available` | `config/modules.tsv` lacks the tool, or Sapelo2 retired it | floating versions; the old counting command changed behaviour silently |
| `measured ... but config/datasets.tsv says ...` | The strand check disagrees with the recorded value | the months when unstranded libraries were counted at `-s 2` |
| `differ in strandedness` / `mix PE and SE` | One study can't share one featureCounts run | fragment and read counts on one scale |
| `Sample IDs repeated within study` | Two libraries share an ID | two columns mapping to one metadata row |
| `found more than once` / `No _2 mate` / `mixes paired-end` | `raw/` is not what it should be | a copy, a partial delivery, a stray file |
| `below 70%` (exit 3) | Assigned/(Assigned+NoFeatures) is low | near 50% everywhere means the strand is wrong |
| `counts are never overwritten` | That study was already counted in this run | silently replacing a matrix someone used |
