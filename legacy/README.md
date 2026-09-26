# Legacy — provenance only, do not run

These scripts produced the published trimmed reads, BAMs and count matrices.
They are kept as the record of what ran. The pipeline in `../bin`, `../slurm`
and `../config` replaces all of them.

## pipeline_2026-08/

The August 2026 pipeline. Its `count_study.sh` produced both published
matrices: `Lauderdale_GRCm38_20260813` (48 libraries) and
`Duncan_GRCm38_20260816` (6 libraries). Its `ADDING_A_BATCH.md` and `PROJECT_LAYOUT.md`
describe that pipeline, not the current one.

## Trimming and alignment scripts

| Script | Produced | Notes |
|---|---|---|
| `trim_2026_Novogene_X202SC26055473.sh` | Trimmed FASTQ, Lauderdale 2026 batch (21 libraries), job 46435579, 25–26 Jun 2026 | Trimmomatic 0.39 (jar path hard-coded); FastQC and MultiQC loaded as "latest available" (MultiQC 1.28 per the report). Post-trim FastQC on R1 only. |
| `trim_Duncan_GSE183742.sh` | Trimmed FASTQ, Duncan GSE183742 (6 libraries) | Identical to the 2026 script apart from `BASE_DIR` and `RAW_DIR`. |
| `align_count_2025_HISAT2_StringTie_featureCounts.sh` | BAMs, Lauderdale 2025 batch (27 libraries), job 38507686, Jul 2025 | Modules hard-coded: HISAT2 2.2.1 (gompi-2022a), SAMtools 1.16.1, StringTie 2.2.1, Subread 2.0.6. Index `genome_snp_tran` (GRCm38), vM25 splice sites. |

## Parameters shared with the current pipeline

- **Trimming.** `PE -phred33 ILLUMINACLIP:<adapters>:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50`. Adapters auto-detected per sample from FastQC overrepresented sequences, with TruSeq3-PE as the fallback.
- **Alignment.** `hisat2 -p 8 --dta --known-splicesite-infile <vM25 splice sites>` against `genome_snp_tran`. That index is transcript-aware but contains no SNPs (see `../config/assemblies.tsv`).

## Superseded, and why

- **featureCounts.** These scripts counted with `-p` but without `--countReadPairs`. Under Subread 2.0.6 that counts reads, not fragments. They also used `-s 2`, but the libraries are unstranded. The published matrices were re-counted from these BAMs by `pipeline_2026-08/count_study.sh`.
- **StringTie and prepDE.** Both ran with the wrong settings, `--rf` on unstranded libraries and prepDE at 75 bp against 150 bp reads. Their output was never used.
- **Output check.** The 2025 script checks `gene_counts.txt` after featureCounts writes `gene_counts_featureCounts.txt`, so it exits before prepDE runs.

## Not included

- **The 2026 alignment script** (`HISAT2_..._PATCH_2.sh`). It used the same HISAT2 parameters; its log confirms HISAT2 2.2.1, SAMtools 1.23.1 and Subread 2.0.6.
- **The 2025 trimming script.** Its MultiQC version, 1.14, matches the fallback in the trimming scripts above, but the script itself has not been located.
- **HISAT2's helper scripts and `prepDE.py3`.** These are distributed with HISAT2 and StringTie under their own licences.

## Tool versions for the 2025 alignment

The versions above come from the script's hard-coded `ml` lines. The script runs under `set -x` and sends stderr to `HISAT2_pipeline_error.log`, so that log, if it survives, records the `ml` commands that actually executed.
