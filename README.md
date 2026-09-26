# pax6-cornea-trigeminal-rnaseq-upstream

FASTQ-to-count-matrix pipeline for the Lauderdale lab's *Pax6* mouse cornea and
trigeminal ganglion RNA-seq, run on the University of Georgia's Sapelo2
cluster. Its count matrices are the shared input to three analyses:

| Repository | Analysis |
|---|---|
| [pax6-cornea-nerve-vascular-transparency](https://github.com/Lauderdale-Lab/pax6-cornea-nerve-vascular-transparency) | Nerve, vascular and transparency study (GEO GSE348661) |
| [pax6-adult-cornea-transcriptome](https://github.com/Lauderdale-Lab/pax6-adult-cornea-transcriptome) | Adult-cornea whole-transcriptome differential expression |
| [pax6-cornea-development-transcriptome](https://github.com/Lauderdale-Lab/pax6-cornea-development-transcriptome) | Developmental differential expression |

James D. Lauderdale, PhD — Department of Cellular Biology, University of
Georgia, Athens, GA 30602, USA.

## In one screen

```bash
bin/pax6 check                                   # config, data, reference, modules
bin/pax6 new-run reprocess_2026-10               # a run = one commit, one genome build
bin/pax6 submit  reprocess_2026-10 Lauderdale_2025 Lauderdale_2026 Duncan_GSE183742
bin/pax6 status  reprocess_2026-10
bin/pax6 count   reprocess_2026-10 Lauderdale    # waits for the strand checks
```

`submit` queues four chained SLURM jobs per dataset:

```
qc_trim (array, one task per library)   FastQC → Trimmomatic → FastQC on both mates
   └─ align (array; task i starts when trim task i succeeds)   HISAT2 → sorted, indexed BAM
        ├─ strand   featureCounts at -s 0/1/2 → measured layout and strandedness
        └─ multiqc  one report per dataset
```

`count` makes one matrix per study in a single featureCounts run, and writes
beside it a `provenance/` folder: the software versions each tool reported,
checksums of the reference files and of every raw FASTQ, a per-library QC
table (reads in, surviving trimming, alignment rate, assignment), the measured
strandedness, and a Methods paragraph filled in from those records. It refuses
to run unless every library has a verified BAM, the strand check measured in
this run agrees with `config/datasets.tsv`, and all datasets in the study share
one strand setting and one layout.

Step-by-step instructions for Sapelo2 are in **[docs/RUNBOOK.md](docs/RUNBOOK.md)**.

## Design

**Data is configuration, not code.** Datasets, genome builds, module versions
and parameters live in `config/`. Adding a dataset is one row in
`config/datasets.tsv`. No script is edited.

**Heavy work on scratch, results on /work.** Runs are processed on
`/scratch`, which is fast and keeps trimmed reads and BAMs off the group
quota. Each count job copies its matrix, provenance and QC to `/work` as soon
as it finishes, so nothing of value depends on scratch surviving.

**A run is one complete, frozen processing.** `new-run` records the pipeline
commit, genome build, module versions and parameters in `runs/<run>/`.
Everything the run produces stays inside that folder, so outputs never mix. If
pipeline code or parameters change, the old run refuses further work and a new
run begins. `datasets.tsv` may still change mid-run, because strand and layout
are recorded as they are measured.

**Every step verifies its product.** Outputs are written under temporary names
and renamed only after they check out: gzip integrity for trimmed reads,
`samtools quickcheck` plus an index for BAMs. A resubmitted job skips what is
already verified and redoes the rest, so a job killed mid-write can't leave a
file that looks finished.

**Nothing is assumed that can be measured.** Layout and strandedness are
measured in every run. The index's contents were checked with `hisat2-inspect`
rather than read off its name.

**Each guard answers a real failure.** They are listed in
[docs/RUNBOOK.md](docs/RUNBOOK.md#why-it-stops).

## Repository layout

```
bin/pax6                  the only command you run
lib/common.sh             shared functions (config tables, runs, discovery, checks)
lib/provenance.sh         the provenance record and Methods draft written with each matrix
slurm/                    one job script per step
  qc_trim.sbatch  align.sbatch  strand.sbatch  multiqc.sbatch  count.sbatch  build_index.sbatch
config/
  datasets.tsv            one row per dataset: study, batch, layout, strand
  assemblies.tsv          genome builds and what their HISAT2 index really contains
  modules.tsv             pinned Lmod modules
  parameters.env          trimming, alignment and counting parameters
  site.env                Sapelo2 paths and SLURM defaults
tools/
  archive_legacy_outputs.sh   move the August 2026 outputs out of the way
  compare_counts.py           compare two count matrices library by library
  make_test_dataset.sh        a few-minute test dataset from real libraries
  verify_md5.sh               check raw FASTQ against vendor checksums
  keep_run.sh                 copy a run's results (optionally BAMs) from scratch to /work
tests/run_tests.sh        75 tests; no cluster needed (stand-ins in tests/stubs/)
docs/                     RUNBOOK.md, LAYOUT.md, PROVENANCE_template.md
legacy/                   the scripts that produced the published data; not for running
```

## Parameters

| Step | Tool (pinned) | Settings |
|---|---|---|
| QC | FastQC 0.12.1 | raw reads and both trimmed mates |
| Trimming | Trimmomatic 0.39 | `ILLUMINACLIP:<adapters>:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50`. Adapters auto-detected from FastQC overrepresented sequences, TruSeq3 otherwise; which was used is recorded per library |
| Alignment | HISAT2 2.2.1, SAMtools 1.18 | `--dta --known-splicesite-infile <GENCODE vM25 splice sites>`; GRCm38 primary assembly; transcript-aware index (splice sites and exons built in, no SNPs) |
| Counting | featureCounts, Subread 2.0.6 | `-p --countReadPairs -B -C -t exon -g gene_id -O -M --primary -Q 10 -s <measured>` |
| Report | MultiQC 1.28 | per dataset |

These trimming, alignment and counting settings are the ones behind the
published matrices, so a reprocessing run is directly comparable with them
(`tools/compare_counts.py`).

## Provenance of the published matrices

| Matrix | Libraries | Produced by |
|---|---|---|
| `Lauderdale_GRCm38_20260813` | 48 (34 cornea, 14 trigeminal) | trimming and alignment by the scripts in `legacy/`; counting by `legacy/pipeline_2026-08/count_study.sh` |
| `Duncan_GRCm38_20260816` | 6 (GSE183742) | as above |

`legacy/README.md` sets out what each old script did, and which of its defects
were corrected before those matrices were made.

## Maintaining it

- **Before every commit:** run `bash tests/run_tests.sh`. GitHub runs the tests
  and shellcheck on every push.
- **Change a module version or parameter:** do it in its own commit and start
  a new run. An existing run keeps what it was created with.
- **New dataset:** add a row to `config/datasets.tsv` and put its FASTQ in
  `data/<dataset>/raw/`. See the runbook.
- **New genome build:** add a row to `config/assemblies.tsv`, place the FASTA
  and matching GENCODE annotation in `reference/<assembly>/`, then run
  `bin/pax6 index <assembly>`.

## Licence

MIT — see `LICENSE`.
