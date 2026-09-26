# /work/jdllab project layout

Last revised: 25 September 2026

```
/work/jdllab/PAX6_RNAseq/
│
├── pipeline/                  scripts + studies.tsv, assemblies.tsv,
│                              module_versions.txt        (version-control this)
├── docs/                      this file, ADDING_A_BATCH.md, decisions, provenance
├── logs/                      SLURM .out/.err, moved here periodically
│
├── reference/                 one subdirectory per genome build
│   ├── GRCm38/                genome FASTA, GENCODE vM25, HISAT2 index
│   └── GRCm39/                genome FASTA, matching GENCODE, HISAT2 index
│
├── data/                      ONE DIRECTORY PER DATASET, IDENTICAL SHAPE
│   ├── Lauderdale_2025/
│   │   ├── PROVENANCE.md      where it came from, when, vendor, accession
│   │   ├── metadata/          sample sheet, vendor QC report, MD5 manifest
│   │   ├── raw/               FASTQ exactly as delivered — never edited
│   │   ├── trimmed/           trimmed FASTQ
│   │   ├── qc/                raw/ trimmed/ adapters/ + multiqc report
│   │   └── align/
│   │       ├── GRCm38/        BAMs + .bai + *_hisat2_summary.txt
│   │       └── GRCm39/        the same libraries against another build
│   ├── Lauderdale_2026/       ← identical shape
│   └── Duncan_GSE183742/      ← identical shape
│
└── counts/                    one matrix per STUDY per ASSEMBLY, dated
    ├── Lauderdale_GRCm38_20260813/
    └── Duncan_GSE183742_GRCm38_20260901/
```

## The four rules

**1. Every dataset has the same internal shape.** `raw/ trimmed/ qc/ align/
metadata/` and a `PROVENANCE.md`. A new dataset is `mkdir` of that shape plus a
row in `studies.tsv` — no decisions, no new directory conventions to remember.
When every dataset looks alike, a script can loop over them and a person can
find anything in one guess.

**2. Datasets are directories; studies are rows in `studies.tsv`.** A dataset
is what arrived from a sequencing centre or a GEO accession. A study is what
you model together. `Lauderdale_2025` and `Lauderdale_2026` are two datasets in
the `Lauderdale` study; `Duncan_GSE183742` is its own study because the allele
differs. Regrouping later means editing one column, not moving files.

**3. Counts live outside the datasets.** A count matrix is built from *all*
libraries of a study at once, so it belongs to no single dataset. The directory
name carries study, assembly and date, so adding a batch produces a new matrix
beside the old rather than replacing it — which is what lets you attribute a
change in results to the thing you changed.

**4. `raw/` is never written to again.** Trimming reads from it and writes
elsewhere. If you ever need to reprocess from scratch, that directory is the
only thing you cannot regenerate. Consider `chmod -R a-w` on it once the
checksums verify.

## Why `align/` is split by assembly

BAMs go to `align/<assembly>/`, so the same libraries can exist against GRCm38
and GRCm39 at once without colliding. `count_study.sh` reads one assembly's
directory only, so a matrix cannot contain a mixture. Selecting a build is
`PAX6_ASSEMBLY=GRCm39` on the submit line; nothing else changes.

## PROVENANCE.md — the file that is always missing

One per dataset, written when the data lands, before anything is run on it:

```markdown
# Lauderdale_2026

Source        : Novogene, project X202SC26055473
Received      : 2026-06-24
Libraries     : 21 (cornea only)
Layout        : 2 x 150 bp paired-end
Strandedness  : unstranded (measured 2026-08-13; the vendor sheet did not say)
Original path : /work/jdllab/Mouse_Pax6_Cornea_2_2026/Novogene_X202SC26055473/01.RawData
Checksums     : metadata/MD5.txt, verified 2026-08-13
Vendor report : metadata/methods_med.pdf
Notes         : Adapter autodetection found nothing in any library; all
                trimmed against TruSeq3-PE.
```

Two of those lines are worth the whole file. **Original path** keeps older
outputs interpretable after a reorganisation — the published count matrix
records BAM paths in its header, and without this note those paths become
mysterious. **Strandedness** records a measured fact that took a full re-count
to establish and that no vendor document states.

## Naming

| Thing | Convention | Example |
|---|---|---|
| Dataset directory | `<Lab>_<year>` or `<Lab>_<accession>` | `Duncan_GSE183742` |
| Study (in `studies.tsv`) | `<Lab>` or the accession | `Lauderdale` |
| Batch (in `studies.tsv`) | year, or the accession | `2026` |
| BAM | `<sample_id>_hisat2.sorted.bam` | `A_CW3_hisat2.sorted.bam` |
| Count matrix dir | `<study>_<assembly>_<date>` | `Lauderdale_GRCm38_20260813` |

No project tags, batch labels or dates inside filenames. That identity lives in
`studies.tsv` and in the `sample_map.tsv` written beside every matrix. The 2026
BAMs carried `_Novogene_X202SC26055473_` in their names and it had to be
stripped with string surgery in R.

Sample IDs must be unique **within a study**. Across studies they may repeat,
since counting never spans studies.

## Space

Raw FASTQ dominates: roughly 10 GB per library at this depth, so about 500 GB
for the 48 existing libraries, plus a similar amount again for BAMs. Check the
group quota before moving anything, and before accepting a new dataset:

```bash
lfs quota -h -g jdllab /work    # or: du -sh /work/jdllab/*
```

Intermediate SAM files are never written — `hisat2` pipes into `samtools sort` —
which is what keeps this tractable. The old pipeline kept a ~200 GB SAM per
library.

## Adding a dataset

```bash
D=/work/jdllab/PAX6_RNAseq/data/NewLab_2027
mkdir -p $D/{raw,trimmed,qc,align,metadata}
cp docs/PROVENANCE_template.md $D/PROVENANCE.md   # then fill it in
# put FASTQ in $D/raw/, verify checksums, then add one row to studies.tsv:
#   NewLab <TAB> 2027 <TAB> auto <TAB> auto <TAB> $D/trimmed <TAB> $D/align
```

Then the usual four jobs: `qc_trim_batch.sh`, `align_batch.sh`,
`strand_check.sh`, `count_study.sh`.

## What this replaces

| Now scattered at | Becomes |
|---|---|
| `/work/jdllab/Mouse_Cornea_and_Trigeminal_May_2025/` | `data/Lauderdale_2025/` |
| `/work/jdllab/Mouse_Pax6_Cornea_2_2026/Novogene_X202SC26055473/01.RawData/` | `data/Lauderdale_2026/raw/` |
| `/work/jdllab/Novogene_X202SC26055473_RNAseq_Map_Counts_Output/RNAseq_Mapping/` | `data/Lauderdale_2026/align/GRCm38/` |
| `/work/jdllab/HISAT2_Mouse_Index/GRCm38/` | `reference/GRCm38/` |
| `/work/jdllab/HISAT2_Mouse_Index/GRCm39_Data/` | `reference/GRCm39/` |
| `/work/jdllab/PAX6_Combined_Counts_20260813/` | `counts/Lauderdale_GRCm38_20260813/` |

Moving within `/work/jdllab` is a rename, not a copy: instant, no extra space,
no risk of a half-copied file. But it *does* invalidate every path recorded in
existing outputs, which is why `PROVENANCE.md` records the original path and
why `studies.tsv` must be updated in the same operation.
