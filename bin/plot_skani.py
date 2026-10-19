#!/usr/bin/env python3
"""
Pairwise ANI and alignment fraction (AF) between accessory plasmids (skani dist).

`split` writes one FASTA per plasmid and a list of them: skani run on one file per
plasmid reproduces the publication exactly, whereas --qi/--ri on a multi-FASTA
gives slightly different values for ~4% of pairs.

`plot` adds AF = min(Align_fraction_ref, Align_fraction_query), i.e. the aligned
fraction of the longer plasmid of each pair, drops self-comparisons, and plots AF vs
ANI as a hexbin density (publication Fig. S2D).

Outputs: plasmid_ani.tsv, ani_vs_af.png / .svg
"""
import argparse
import os
import sys


def split(args):
    os.makedirs("plasmids", exist_ok=True)
    names, out = [], None
    with open(args.fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                if name in names:
                    sys.exit(f"ERROR: duplicated sequence {name}")
                names.append(name)
                if out:
                    out.close()
                out = open(os.path.join("plasmids", f"{name}.fasta"), "w")
                out.write(f">{name}\n")
            elif out:
                out.write(line)
    if out:
        out.close()
    with open("plasmids.list", "w") as fh:
        fh.writelines(os.path.abspath(os.path.join("plasmids", f"{n}.fasta")) + "\n" for n in names)
    print(f"{len(names)} plasmids")


def plot(args):
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.read_csv(args.skani, sep="\t")
    df["Ref_name"] = df["Ref_name"].str.split().str[0]
    df["Query_name"] = df["Query_name"].str.split().str[0]
    df = df[df["Ref_name"] != df["Query_name"]].copy()
    df["AF"] = df[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
    cols = ["Ref_name", "Query_name", "ANI", "Align_fraction_ref", "Align_fraction_query", "AF"]
    df[cols].to_csv("plasmid_ani.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7, 6))
    hb = ax.hexbin(df["ANI"], df["AF"], gridsize=60, cmap="viridis", mincnt=1, bins="log")
    fig.colorbar(hb, ax=ax).set_label("Plasmid pairs (log scale)")
    ax.set_xlabel("ANI (%)")
    ax.set_ylabel("Alignment fraction of the longer plasmid (%)")
    ax.set_title(f"AF vs ANI, accessory plasmids all-vs-all (n = {len(df)} pairs)")
    ax.set_ylim(0, 100)
    if len(df):
        ax.set_xlim(df["ANI"].min() - 1, 100)
    fig.savefig("ani_vs_af.png", dpi=300, bbox_inches="tight")
    fig.savefig("ani_vs_af.svg", bbox_inches="tight")
    print(f"{len(df)} plasmid pairs with ANI reported by skani")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("split")
    a.add_argument("--fasta", required=True)
    b = sub.add_parser("plot")
    b.add_argument("--skani", required=True, help="skani dist output")
    args = p.parse_args()
    split(args) if args.cmd == "split" else plot(args)


if __name__ == "__main__":
    main()
