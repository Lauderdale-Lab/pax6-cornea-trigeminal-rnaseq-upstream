# PAX6 mouse cornea RNA-seq — upstream processing (FASTQ → gene counts)

Shell pipeline that takes raw paired-end FASTQ files to gene-level count
matrices for the Lauderdale lab's *Pax6*<sup>Sey-Neu</sup>/+ mouse cornea and
trigeminal ganglion RNA-seq, and for the independent dataset GEO GSE183742
used for cross-dataset replication.

The same count matrices are analysed in three papers, each with its own
analysis repository:

| Paper | Analysis repository |
|---|---|
| Nerve and vascular abnormalities precede loss of corneal transparency in *Pax6*-haploinsufficient mice (Mohan & Lauderdale) | [pax6-cornea-nerve-vascular-transparency](https://github.com/Lauderdale-Lab/pax6-cornea-nerve-vascular-transparency) |
| Adult cornea whole-transcriptome analysis (in preparation) | [pax6-adult-cornea-transcriptome](https://github.com/Lauderdale-Lab/pax6-adult-cornea-transcriptome) |
| Developmental cornea transcriptome analysis (in preparation) | [pax6-cornea-development-transcriptome](https://github.com/Lauderdale-Lab/pax6-cornea-development-transcriptome) |

Repository: <https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream> ·
Archived release: `<ZENODO_DOI>` · Contact: James D. Lauderdale, <jdlauder@uga.edu>

Licensed under the MIT License — see [`LICENSE`](LICENSE).

---

## What it does

1. **QC and trimming** — FastQC on raw reads; Trimmomatic 0.39 (PE, `-phred33`,
   `ILLUMINACLIP:<adapters>:2:30:10:2:True LEADING:3 TRAILING:3
   SLIDINGWINDOW:4:20 MINLEN:50`); FastQC again on **both** trimmed mates;
   MultiQC.
2. **Alignment** — HISAT2 2.2.1 to GRCm38 (`genome_snp_tran` index), with
   GENCODE vM25 splice sites supplied at alignment time
   (`--known-splicesite-infile`), piped straight into `samtools sort`.
3. **Strandedness is measured, not assumed** — `strand_check.sh` counts a
   subset of BAMs at `-s 0`, `-s 1` and `-s 2` and applies a written rule.
4. **Counting** — featureCounts (Subread 2.0.6) against GENCODE vM25, run once
   per study over every library in a single invocation, counting fragments
   (`-p --countReadPairs -B -C`).

Transcript assembly (StringTie) and `prepDE.py` are not part of this pipeline.

## Design: datasets are data, not code

Adding a sequencing batch or a new study means adding a row to a manifest, not
editing a script.

| File | Holds |
|---|---|
| `pipeline/studies.tsv` | one row per batch: study, batch, layout, strand, trimmed and mapping directories |
| `pipeline/assemblies.tsv` | one row per genome build: reference directory, FASTA, annotation, index |
| `pipeline/module_versions.txt` | the pinned Lmod module for each tool |

A **dataset** is files that arrived together; a **study** is libraries modelled
together. The two Lauderdale batches are two datasets and one study;
GSE183742 is one of each. Counting refuses to mix studies, genome builds,
library layouts (paired vs single end) or strandedness in one matrix.

See `docs/ADDING_A_BATCH.md` for the step-by-step runbook and
`docs/DIRECTORY_LAYOUT.md` for the directory structure.

## Datasets processed

| Dataset | Libraries | Layout | Strandedness (measured) | Count matrix |
|---|---|---|---|---|
| Lauderdale_2025 | 27 (13 cornea, 14 trigeminal ganglion) | PE 150 | unstranded (`-s 0`) | Lauderdale, 48 libraries |
| Lauderdale_2026 (Novogene X202SC26055473) | 21 cornea | PE 150 | unstranded (`-s 0`) | Lauderdale, 48 libraries |
| Duncan_GSE183742 | 6 adult cornea | PE 101 | reverse (`-s 2`) | Duncan, 6 libraries |

Raw reads for the Lauderdale libraries are deposited at GEO `<ACCESSION>`.
GSE183742 (Krishnan, Faranda, Novo, Wang & Duncan) was reprocessed from SRA
SRP336260 with the same reference and annotation.

The Lauderdale count matrix used by the analysis repositories is
`gene_counts_featureCounts_all48.txt`
(MD5 `67735b701252702ff4d1be50366e30e6`).

## Provenance: which code produced the published files

- **BAMs.** The Lauderdale BAMs were produced before this pipeline was
  written, by the scripts in [`legacy/`](legacy/). They are kept for
  provenance only and must not be run. The alignment reference and parameters
  are identical to those in `align_batch.sh`. HISAT2, SAMtools and Subread
  versions for the 2026 batch are recorded in its job logs. The 2025 batch's
  logs did not record tool versions.
- **Counts.** `count_study.sh` reproduces the published Lauderdale matrix: all
  48 libraries, with every library's assigned fraction identical to the
  original run.

## Findings that shaped the pipeline

- **The Lauderdale libraries are not strand-specific.** Counting at `-s 1` and
  `-s 2` assigns the same ~35% of fragments; `-s 0` assigns 65–78%. This holds
  in both batches. GSE183742, by contrast, is reverse-stranded.
- **`-p` alone does not count fragments in Subread ≥ 2.0.2.**
  `--countReadPairs` is required; without it featureCounts silently assigns
  in single-end mode and `-B`/`-C` have no effect.
- **GSE183742 aligns at 51–55%.** The unaligned reads are dominated by 45S
  pre-rRNA spacer sequence (ribo-depleted library), which falls in the rDNA
  array that is collapsed in the primary assembly. Assigned fractions among
  aligned reads are normal.

## Running it

Written for a SLURM cluster with Lmod modules (UGA Sapelo2). On another
system, edit the three manifest files; no script should need changing.

```bash
cd <data root>/pipeline
sbatch qc_trim_batch.sh   <study> <batch>
sbatch align_batch.sh     <study> <batch>
sbatch strand_check.sh    <study> <batch>
# record the measured strandedness in studies.tsv, then:
sbatch count_study.sh     <study>
```

The genome build is chosen with `PAX6_ASSEMBLY` (default: the first row of
`assemblies.tsv`, GRCm38). Every published result uses GRCm38 with GENCODE
vM25.

## Requirements

Bash, SLURM, Lmod, and the modules pinned in `module_versions.txt`: HISAT2
2.2.1, SAMtools 1.23.1, Subread 2.0.6, plus Trimmomatic 0.39, FastQC and
MultiQC.

## Citing

See `CITATION.cff`, or use the "Cite this repository" button on GitHub.

## Licence

MIT — see [`LICENSE`](LICENSE).
