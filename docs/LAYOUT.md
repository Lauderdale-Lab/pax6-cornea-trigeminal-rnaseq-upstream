# Layout on /work

```
/work/jdllab/PAX6_RNAseq/
├── pipeline/                  this repository (git clone); the only code
├── reference/<assembly>/      genome FASTA, GENCODE annotation, HISAT2 index
├── data/<dataset>/            INPUTS ONLY, one folder per delivery
│   ├── PROVENANCE.md          source, dates, measured strandedness, original path
│   ├── metadata/              vendor reports, MD5 manifests, optional samples.tsv
│   └── raw/                   FASTQ exactly as delivered; read-only
├── runs/<run>/                RESULTS KEPT: copied here from scratch by every count
│                              job (everything below except trimmed/ and BAMs;
│                              BAMs too after tools/keep_run.sh --with-bams)
│   ├── RUN_INFO.tsv           commit, genome build, creator, date
│   ├── modules.tsv            module versions frozen for this run
│   ├── parameters.env         parameters frozen for this run
│   ├── jobs.tsv               every job submitted, by dataset
│   ├── logs/                  SLURM output of every job
│   ├── reference/             splice sites extracted for this run
│   ├── provenance/versions/   what each job's tools reported, one file per job
│   ├── <dataset>/
│   │   ├── samples.tsv        the libraries, in array-task order
│   │   ├── checksums/         MD5 of each library's raw FASTQ
│   │   ├── trimmed/           <id>_1.fq.gz, <id>_2.fq.gz, <id>.done
│   │   ├── qc/                raw/ trimmed/ FastQC; adapters/; trimmomatic/ logs
│   │   ├── align/             <id>.bam, .bai, <id>.hisat2.txt
│   │   ├── strand/            strand_table.tsv, RESULT.tsv
│   │   └── multiqc_<dataset>.html
│   └── counts/<study>/
│       ├── gene_counts.tsv    gene_id x sample_id: the matrix analyses read
│       ├── featureCounts.txt  raw featureCounts output (+ .summary)
│       ├── sample_map.tsv     sample_id, study, dataset, batch, BAM
│       ├── assignment_qc.tsv  per-library assignment rates
│       ├── PROVENANCE.tsv     run, commit, build, strand, command, matrix MD5
│       └── provenance/        software_versions, reference (MD5s), inputs (raw
│                              FASTQ MD5s), library_qc, strand, METHODS.md
└── archive/                   superseded outputs, never read by the pipeline

/scratch/$USER/PAX6_RNAseq/
└── runs/<run>/                PROCESSING: the full run, same structure as above,
                               including trimmed/ and BAMs. Purged by Sapelo2,
                               so nothing here is permanent.
```

## Rules

1. **Inputs and outputs never share a folder.** `data/` is written once, when
   a delivery arrives. Everything computed goes to `runs/`.
2. **Every dataset has the same shape.** `raw/`, `metadata/` and a
   `PROVENANCE.md`, so scripts and people find things in one guess.
3. **A dataset is a folder; a study is a column.** Datasets are what arrived
   together. Studies are what gets modelled together, set in
   `config/datasets.tsv`. Regrouping means editing a column, not moving files.
4. **Runs are never edited in place.** A new build, module version or
   parameter set means a new run beside the old one, which is what lets a
   change in results be traced to the change that caused it.
5. **Space.** Raw FASTQ is about 10 GB per library at this depth, and a run
   adds about 15 GB per library of trimmed reads and BAMs. That is why runs
   are processed on scratch. On /work, a kept run without BAMs is a few GB;
   with BAMs, add about 3 GB per library. Check the quota with
   `lfs quota -h -g jdllab /work`.

## Naming

| Thing | Convention | Example |
|---|---|---|
| Dataset | `<Lab>_<year>` or `<Lab>_<accession>` | `Lauderdale_2026`, `Duncan_GSE183742` |
| Study | `<Lab>` | `Lauderdale` |
| Run | purpose and date | `reprocess_2026-10`, `test_20261001` |
| Sample ID | from the file name, or `metadata/samples.tsv` | `A_CW3` |
