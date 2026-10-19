#!/usr/bin/env python3
"""
Trinucleotide (3-mer) composition of accessory plasmids vs the main replicons.

`prepare` pools the main replicons (contigs tagged Chromosome / pSymA / pSymB) and
the accessory plasmids of the samples (PLSDB plasmids excluded) into one FASTA, and
writes a table of the same sequences in the same order (Kmertools output has no IDs).

`analyse` reads the Kmertools frequencies (`kmertools comp oligo -k 3`: normalised
frequencies of the 32 canonical 3-mers, one row per sequence in FASTA order),
computes Bray-Curtis distances and a UMAP (precomputed distances, random_state 42,
default n_neighbors 15 and min_dist 0.1). The publication notebook passed
n_neighbors=10 to fit_transform(), where UMAP ignores it, so the defaults were used
there too. Plots:
  kmer_umap          UMAP of 3-mer composition coloured by replicon type
  replicon_length    length per replicon type
  replicon_gc        GC content per replicon type
Length and GC of accessory plasmids are compared with each main replicon type
(two-sided Mann-Whitney U, Bonferroni-corrected).

Outputs (analyse): kmer_3mer_frequencies.tsv, kmer_umap.tsv, kmer_stats.tsv, figures (PNG + SVG)
"""
import argparse
import csv
import os
import sys

ACCESSORY = "Accessory plasmid"
# Kmertools 1.0.0 column order for k = 3 (canonical 3-mers)
KMERS3 = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "ATA", "ATC", "ATG", "CAA",
          "CAC", "CAG", "CCA", "CCC", "CCG", "CGA", "CGC", "CTA", "CTC", "GAA", "GAC", "GCA", "GCC", "GGA", "GTA",
          "TAA", "TCA"]
ORDER = ["Chromosome", "pSymA", "pSymB", ACCESSORY]


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split(None, 1)[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def prepare(args):
    tags = {}
    for path in args.tags:
        with open(path) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                tags[f[0]] = f[1]
    exclude = set()
    if args.exclude:
        with open(args.exclude) as fh:
            exclude = {l.strip() for l in fh if l.strip()}

    seen = set()
    with open("kmer_input.fasta", "w") as fa, open("kmer_input.tsv", "w") as tab:
        tab.write("sequence\ttype\tlength\tgc\n")

        def emit(name, seq, kind):
            if name in seen:
                sys.exit(f"ERROR: duplicated sequence {name}")
            seen.add(name)
            s = seq.upper()
            gc = 100 * (s.count("G") + s.count("C")) / len(s) if s else 0
            fa.write(f">{name}\n{seq}\n")
            tab.write(f"{name}\t{kind}\t{len(seq)}\t{gc:.2f}\n")

        for path in sorted(args.main):
            for name, seq in read_fasta(path):
                emit(name, seq, tags.get(name, "Other"))
        n_pls = 0
        for name, seq in read_fasta(args.plasmids):
            if name not in exclude:
                emit(name, seq, ACCESSORY)
                n_pls += 1
    print(f"{len(seen) - n_pls} main replicons + {n_pls} accessory plasmids")


def analyse(args):
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import umap
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import mannwhitneyu

    meta = pd.read_csv(args.table, sep="\t", dtype={"sequence": str})
    freqs = np.loadtxt(args.kmers, ndmin=2)
    if freqs.shape != (len(meta), len(KMERS3)):
        sys.exit(f"ERROR: Kmertools gave {freqs.shape} values for {len(meta)} sequences x {len(KMERS3)} 3-mers")
    table = pd.concat([meta, pd.DataFrame(freqs, columns=KMERS3)], axis=1)
    table.to_csv("kmer_3mer_frequencies.tsv", sep="\t", index=False, float_format="%.6f")

    dist = squareform(pdist(freqs, metric="braycurtis"))
    n_neighbors = min(15, len(meta) - 1)  # small test sets have fewer than 15 sequences
    emb = umap.UMAP(random_state=42, metric="precomputed", n_neighbors=n_neighbors).fit_transform(dist)
    meta["UMAP1"], meta["UMAP2"] = emb[:, 0], emb[:, 1]
    meta.to_csv("kmer_umap.tsv", sep="\t", index=False)

    order = [o for o in ORDER if o in set(meta["type"])]
    palette = dict(zip(order, sns.color_palette("colorblind", len(order))))
    sns.set_theme(style="ticks")

    fig, ax = plt.subplots(figsize=(6, 5))
    for kind in order:
        sub = meta[meta["type"] == kind]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=12, alpha=0.85, color=palette[kind], label=f"{kind} (n = {len(sub)})",
                   linewidths=0)
    ax.set_title("UMAP of 3-mer composition (Bray-Curtis)")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig("kmer_umap.png", dpi=300, bbox_inches="tight")
    fig.savefig("kmer_umap.svg", bbox_inches="tight")
    plt.close(fig)

    # Accessory plasmids vs each main replicon type (two-sided Mann-Whitney, Bonferroni)
    stats = []
    mains = [o for o in order if o != ACCESSORY]
    for col in ("length", "gc"):
        acc = meta.loc[meta["type"] == ACCESSORY, col]
        for kind in mains:
            other = meta.loc[meta["type"] == kind, col]
            if len(acc) and len(other):
                u, p = mannwhitneyu(acc, other, alternative="two-sided")
                stats.append({"variable": col, "group1": ACCESSORY, "group2": kind, "n1": len(acc), "n2": len(other),
                              "U": u, "p": p, "p_bonferroni": min(1.0, p * len(mains))})
    pd.DataFrame(stats).to_csv("kmer_stats.tsv", sep="\t", index=False)

    for col, label, name in (("length", "Length (bp)", "replicon_length"), ("gc", "GC (%)", "replicon_gc")):
        fig, ax = plt.subplots(figsize=(5, 4.5))
        sns.violinplot(data=meta, x="type", y=col, hue="type", order=order, hue_order=order, palette=palette,
                       cut=0, legend=False, ax=ax)
        title = label.split(" (")[0] + " by replicon type"
        s = [r for r in stats if r["variable"] == col]
        if s:
            title += "\n" + "; ".join(f"vs {r['group2']}: p = {r['p_bonferroni']:.2g}" for r in s)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel(label)
        ax.tick_params(axis="x", rotation=20)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{name}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{name}.svg", bbox_inches="tight")
        plt.close(fig)
    print(f"{len(meta)} sequences: " + ", ".join(f"{k} {int((meta['type'] == k).sum())}" for k in order))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("prepare")
    a.add_argument("--main", nargs="+", required=True, help="<sample>.main_replicons.fasta files")
    a.add_argument("--tags", nargs="+", required=True, help="<sample>.replicon_tags.tsv files")
    a.add_argument("--plasmids", required=True, help="accessory plasmids FASTA")
    a.add_argument("--exclude", help="IDs to leave out (e.g. PLSDB plasmids)")
    b = sub.add_parser("analyse")
    b.add_argument("--kmers", required=True, help="kmertools comp oligo output")
    b.add_argument("--table", required=True, help="kmer_input.tsv from prepare")
    args = p.parse_args()
    prepare(args) if args.cmd == "prepare" else analyse(args)


if __name__ == "__main__":
    main()
