# Reprocessing check, September 2026

All 54 libraries were reprocessed from raw FASTQ with this pipeline (run
`reprocess_2026-09`, commit `9182fd5`, 26 September 2026). The resulting
count matrices were compared with the matrices analyzed in the published
work, using `tools/compare_counts.py`.

| Study | Libraries | Result |
|---|---|---|
| Duncan (GSE183742) | 6 | identical, every gene in every library |
| Lauderdale (GSE348661) | 48 | library totals within 31 fragments (of 29–59 million); at most 26 fragments for any gene; Pearson r of log2(count + 1) > 0.99999 in every library |

The run used the same trimming, alignment and counting parameters, reference,
index and program versions as the original processing. The only version
differences were SAMtools 1.18 (sorting and indexing only) and FastQC 0.12.1.

Files, per study:

- `compare_with_published.tsv`: the per-library comparison.
- `provenance/`: the record written by the count job:
  - software versions as reported;
  - reference and raw-FASTQ checksums;
  - per-library QC;
  - measured strandedness;
  - the generated Methods draft.

The Duncan `software_versions.tsv` has a second FastQC line holding a JVM
warning, not a version. Only FastQC 0.12.1 ran. The pipeline has since been
fixed to ignore such lines.
