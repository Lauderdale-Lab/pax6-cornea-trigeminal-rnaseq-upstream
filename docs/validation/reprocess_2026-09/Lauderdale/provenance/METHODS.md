# Methods draft — Lauderdale, run reprocess_2026-09

Generated 2026-09-26T12:19:27-04:00 from the recorded provenance of this matrix. Every
number and version below is taken from files in this folder; check the
wording, not the facts.

## RNA-seq processing

Raw reads from 48 libraries were assessed with FastQC 0.12.1 before and after
trimming (both mates). Adapters and low-quality bases were removed with
Trimmomatic 0.39 (ILLUMINACLIP:<adapters>:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50).
Adapter sequences were taken from each library's FastQC overrepresented
sequences where present (0 of 48 libraries) and otherwise from the
Trimmomatic TruSeq3 set. 96.7–98.5% of reads or read pairs
survived trimming (median 97.2%).

Trimmed reads were aligned to the GRCm38 primary assembly
(GRCm38.primary_assembly.genome.fa) with HISAT2 2.2.1 (--dta), using
a transcript-aware HISAT2 index (splice sites and exons from the annotation built in; no SNPs) and supplying 284754
known splice junctions extracted from gencode.vM25.annotation.gtf with
--known-splicesite-infile. Alignments were sorted and indexed with
SAMtools 1.18. Overall alignment rates were 86.7–97.8% (median 91.0%).

Library strandedness was measured, not assumed, by counting a subset of
libraries at each featureCounts strand setting (fragments assigned at -s 0 / 1 / 2 — Lauderdale_2025: 65.5% / 35.1% / 34.2%; Lauderdale_2026: 73.4% / 37.8% / 38.3%). All libraries in
this study were unstranded (featureCounts -s 0).

Fragments were counted per gene with featureCounts (Subread 2.0.6)
(-s 0 -p --countReadPairs -B -C -t exon -g gene_id -O -M --primary -Q 10) against
gencode.vM25.annotation.gtf, all 48 libraries in a single invocation.
65.1–78.4% of fragments were assigned to genes (median 69.2%);
Assigned/(Assigned + Unassigned_NoFeatures) was 81.1–96.9% (median 94.3%).
QC reports were compiled with MultiQC 1.28.

Processing used pipeline commit 9182fd5eb291
(https://github.com/Lauderdale-Lab/pax6-cornea-trigeminal-rnaseq-upstream).
Raw-file checksums: provenance/inputs.tsv. Reference checksums:
provenance/reference.tsv. Per-library values: provenance/library_qc.tsv.
