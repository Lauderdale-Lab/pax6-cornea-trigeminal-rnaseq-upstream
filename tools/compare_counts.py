#!/usr/bin/env python3
"""Compare two gene-count matrices library by library.

    tools/compare_counts.py NEW OLD [--out report.tsv]

Built to check a reprocessing run against the published matrices, e.g.

    tools/compare_counts.py \\
        $PAX6_ROOT/runs/<run>/counts/Lauderdale/gene_counts.tsv \\
        $PAX6_ROOT/archive/2026-08_pipeline/counts/Lauderdale_GRCm38_20260813/gene_counts_featureCounts.txt

Either file may be raw featureCounts output (a '# Program' line, then Geneid,
Chr, Start, End, Strand, Length and one column per BAM path) or a tidy matrix
(gene_id then one column per sample). Column names are reduced to sample IDs
by dropping the directory, '.bam', '_hisat2.sorted' and the old
'_Novogene_<project>' tag, so both naming schemes line up.

Per library it reports total counts in each matrix, their ratio, the fraction
of genes with identical counts, the largest absolute difference, and the
Pearson correlation of log2(count + 1). Standard library only: runs with any
python3.
"""
import argparse
import math
import re
import sys


def sample_id(name):
    s = name.rsplit("/", 1)[-1]
    s = re.sub(r"\.bam$", "", s)
    s = re.sub(r"_hisat2\.sorted$", "", s)
    s = re.sub(r"_hisat2$", "", s)
    s = re.sub(r"_Novogene_.*$", "", s)
    return s


def read_matrix(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#")]
    header = lines[0].split("\t")
    first = 6 if header[:6] == ["Geneid", "Chr", "Start", "End", "Strand", "Length"] else 1
    samples = [sample_id(h) for h in header[first:]]
    if len(set(samples)) != len(samples):
        sys.exit(f"{path}: sample IDs repeat after name clean-up: {samples}")
    counts = {s: {} for s in samples}
    for line in lines[1:]:
        f = line.split("\t")
        gene = f[0].split(".")[0]          # drop the version suffix (ENSMUSG...18)
        for s, v in zip(samples, f[first:]):
            counts[s][gene] = counts[s].get(gene, 0) + int(float(v))
    return counts


def pearson(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("new")
    ap.add_argument("old")
    ap.add_argument("--out", help="also write the per-library table here")
    a = ap.parse_args()

    new, old = read_matrix(a.new), read_matrix(a.old)
    both = sorted(set(new) & set(old))
    only_new, only_old = sorted(set(new) - set(old)), sorted(set(old) - set(new))
    print(f"libraries: {len(both)} in both, {len(only_new)} only in NEW, {len(only_old)} only in OLD")
    if only_new:
        print("  only in NEW:", " ".join(only_new))
    if only_old:
        print("  only in OLD:", " ".join(only_old))
    if not both:
        sys.exit("No libraries in common -- check the sample naming.")

    genes_new, genes_old = set(new[both[0]]), set(old[both[0]])
    genes = sorted(genes_new & genes_old)
    print(f"genes: {len(genes)} in both, {len(genes_new - genes_old)} only in NEW, {len(genes_old - genes_new)} only in OLD")

    rows = []
    for s in both:
        x = [new[s][g] for g in genes]
        y = [old[s][g] for g in genes]
        tn, to = sum(x), sum(y)
        same = sum(1 for p, q in zip(x, y) if p == q) / len(genes)
        maxd, gmax = max(((abs(p - q), g) for p, q, g in zip(x, y, genes)), default=(0, "-"))
        r = pearson([math.log2(v + 1) for v in x], [math.log2(v + 1) for v in y])
        rows.append((s, tn, to, tn / to if to else float("nan"), same, maxd, gmax, r))

    hdr = ("sample_id", "total_new", "total_old", "ratio_new_old", "identical_genes", "max_abs_diff", "gene_max_diff", "pearson_log2")
    fmt = "{:<14} {:>12} {:>12} {:>9} {:>9} {:>12} {:<22} {:>9}"
    print()
    print(fmt.format(*hdr))
    for s, tn, to, ratio, same, maxd, gmax, r in rows:
        print(fmt.format(s, tn, to, f"{ratio:.4f}", f"{same:.1%}", maxd, gmax, f"{r:.5f}"))

    exact = [row[0] for row in rows if row[4] == 1.0]
    print()
    print(f"{len(exact)} of {len(rows)} libraries identical gene for gene.")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\t".join(hdr) + "\n")
            for row in rows:
                fh.write("\t".join(str(v) for v in row) + "\n")
        print(f"table written to {a.out}")


if __name__ == "__main__":
    main()
