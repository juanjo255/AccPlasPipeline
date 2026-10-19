#!/usr/bin/env python3
"""
Plasmid copy number from long-read depth (publication method).

Copy number = median per-base depth of the plasmid / median per-base depth of the
chromosome of the same assembly. The first and last --end-trim bp of each contig
are ignored (reads cannot span the junction of a circular contig) when the contig is
longer than 4 x --end-trim. The chromosome is the longest contig tagged Chromosome
(else the longest contig). pSymA and pSymB are reported too as a sanity check
(expected ~1).

Input depth: `samtools depth -a` output (contig, position, depth), file or '-' for
stdin, possibly gzip-compressed. Contig IDs are renamed as in the contig length
filter (`_<sample>` added if missing).

Output: <sample>.copy_number.tsv
"""
import argparse
import gzip
import sys

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--depth", required=True, help="samtools depth -a output ('-' = stdin; .gz accepted)")
    p.add_argument("--tags", required=True, help="<sample>.replicon_tags.tsv")
    p.add_argument("--plasmids", required=True, help="IDs of accessory plasmids (one per line)")
    p.add_argument("--sample", required=True)
    p.add_argument("--add-suffix", choices=["true", "false"], required=True)
    p.add_argument("--end-trim", type=int, default=1000)
    return p.parse_args()


def depth_stats(depth, end_trim):
    n = len(depth)
    core = depth[end_trim:n - end_trim] if n > 4 * end_trim else depth
    q1, med, q3 = np.percentile(core, [25, 50, 75])
    return {"length": n, "median_depth": med, "mean_depth": core.mean(), "q1_depth": q1, "q3_depth": q3,
            "frac_zero_depth": float((depth == 0).mean())}


def main():
    args = parse_args()
    suffix = f"_{args.sample}"

    def rename(c):
        return c if args.add_suffix == "false" or c.endswith(suffix) else c + suffix

    with open(args.tags) as fh:
        next(fh)
        tags = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tags[f[0]] = (f[1], int(f[2]))
    with open(args.plasmids) as fh:
        plasmids = {l.strip() for l in fh if l.strip()}
    plasmids = {p for p in plasmids if p in tags}  # this sample's plasmids

    chrom = [c for c, (t, _) in tags.items() if t == "Chromosome"]
    chrom_source = "tag" if chrom else "longest_contig"
    chrom_contig = max(chrom or tags, key=lambda c: tags[c][1])
    wanted = plasmids | {c for c, (t, _) in tags.items() if t in ("pSymA", "pSymB")} | {chrom_contig}

    # Collect depths only for the contigs needed (depth tables are large)
    opener = gzip.open if args.depth.endswith(".gz") else open
    fh = sys.stdin if args.depth == "-" else opener(args.depth, "rt")
    values = {c: [] for c in wanted}
    for line in fh:
        c, _, d = line.rstrip("\n").split("\t")[:3]
        c = rename(c)
        if c in values:
            values[c].append(int(d))
    depth = {c: np.asarray(v, dtype=np.int32) for c, v in values.items() if v}
    if chrom_contig not in depth:
        sys.exit(f"ERROR: [{args.sample}] chromosome {chrom_contig} not in the depth table (BAM/assembly mismatch?)")
    for c, arr in depth.items():
        if len(arr) != tags[c][1]:
            sys.exit(f"ERROR: [{args.sample}] {c} has {len(arr)} depth positions but is {tags[c][1]} bp "
                     "(BAM not made against this assembly?)")
    chrom_stats = depth_stats(depth[chrom_contig], args.end_trim)

    cols = ["plasmid", "sample", "type", "length", "median_depth", "mean_depth", "q1_depth", "q3_depth",
            "frac_zero_depth", "chrom_contig", "chrom_source", "chrom_median_depth", "chrom_mean_depth",
            "copy_number", "copy_number_mean", "note"]
    with open(f"{args.sample}.copy_number.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        rows = sorted((("accessory" if c in plasmids else tags[c][0]), c) for c in wanted if c != chrom_contig or c in plasmids)
        for kind, c in rows:
            if c not in depth:
                out.write("\t".join([c, args.sample, kind] + [""] * 12 + ["contig_not_in_depth"]) + "\n")
                continue
            s = depth_stats(depth[c], args.end_trim)
            cn = s["median_depth"] / chrom_stats["median_depth"] if chrom_stats["median_depth"] > 0 else float("nan")
            cnm = s["mean_depth"] / chrom_stats["mean_depth"] if chrom_stats["mean_depth"] > 0 else float("nan")
            vals = [c, args.sample, kind, s["length"], s["median_depth"], s["mean_depth"], s["q1_depth"], s["q3_depth"],
                    s["frac_zero_depth"], chrom_contig, chrom_source, chrom_stats["median_depth"],
                    chrom_stats["mean_depth"], cn, cnm, ""]
            out.write("\t".join(str(v) for v in vals) + "\n")
    print(f"[{args.sample}] chromosome {chrom_contig} median depth {chrom_stats['median_depth']:g}; "
          f"{len(plasmids)} accessory plasmids")


if __name__ == "__main__":
    main()
