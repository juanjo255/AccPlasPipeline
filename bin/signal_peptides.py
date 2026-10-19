#!/usr/bin/env python3
"""
Signal peptide-bearing protein families on the accessory plasmids (publication Fig. 9D).

Port of amp_prediciton.ipynb. SignalP 4.1 (gram-, -s notm) summary rows (all 'SIGNAL'
hits) are linked to the proteins (<locus>__<plasmid> FASTA headers, Bakta product as
description) and to the Panaroo gene families of the plasmid-only pangenome (locus tag
-> family through the gene_presence_absence.csv cells).

  - Plasmid / cluster coverage (paper: "93% of plasmid clusters, 395 plasmids"): plasmids
    with at least one signal-peptide protein in a family whose product is not
    conjugation-related, as in the notebook.
  - Family ranking (Fig. 9D): CAZymes (dbCAN, kept by >= min tools) are removed first,
    families are ranked by the number of plasmids carrying them, conjugation-related
    products ('conjuga') are dropped, and each family is labelled with the product of its
    longest member.

Outputs: signal_peptide_proteins.tsv, signal_peptide_families.tsv,
         signal_peptide_summary.tsv, signal_peptide_families.png/.svg
"""
import argparse
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--signalp", required=True, help="SignalP 4.1 summary (-f summary, GFF-like)")
    p.add_argument("--proteins", required=True, help="plasmid proteins FASTA (<locus>__<plasmid> product)")
    p.add_argument("--gpa", required=True, help="Panaroo gene_presence_absence.csv of the accessory plasmids")
    p.add_argument("--plasmids", required=True, help="plasmid_summary.tsv (sequence, cluster)")
    p.add_argument("--cazymes", help="dbCAN CAZymes kept by the tool filter (cazymes.tsv protein_id, or overview 'Gene ID')")
    p.add_argument("--top", type=int, default=15, help="families shown in the bar plot")
    return p.parse_args()


def main():
    args = parse_args()
    import pandas as pd

    # proteins: name, description, length
    names, descr, lens = [], [], []
    with open(args.proteins) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].rstrip("\n").split(" ", 1)
                names.append(f[0])
                descr.append(f[1] if len(f) > 1 else "")
                lens.append(0)
            elif names:
                lens[-1] += len(line.strip())
    prot = pd.DataFrame({"protein": names, "product": descr, "length": lens})
    prot["locus"] = prot["protein"].str.split("__").str[0]
    prot["plasmid"] = prot["protein"].str.split("__").str[1]

    sp = pd.read_csv(args.signalp, sep="\t", comment="#", header=None,
                     names=["protein", "method", "feature", "start", "end", "score", "na1", "na2", "pred"])
    sp = sp[sp["feature"] == "SIGNAL"]
    df = sp[["protein", "end", "score"]].rename(columns={"end": "cleavage_after"}).merge(prot, on="protein", how="left")

    # Panaroo locus tag -> family
    fam = {}
    with open(args.gpa) as fh:
        r = csv.reader(fh)
        next(r)
        for row in r:
            for cell in row[3:]:
                for token in cell.split(";"):
                    token = token.strip()
                    if token:
                        fam[token] = row[0]
    df["family"] = df["locus"].map(fam)
    df.sort_values("protein").to_csv("signal_peptide_proteins.tsv", sep="\t", index=False)

    conj = df["product"].str.contains("conjuga", case=False, na=False)
    cluster_of = pd.read_csv(args.plasmids, sep="\t", dtype=str).set_index("sequence")["cluster"]
    covered = df[df["family"].notna() & ~conj]["plasmid"].unique()
    n_cl_total = cluster_of.nunique()
    n_cl = cluster_of.reindex(covered).dropna().nunique()

    # Fig. 9D: CAZymes removed, families ranked by number of plasmids
    fam_df = df[df["family"].notna()]
    n_caz = 0
    if args.cazymes:
        caz_df = pd.read_csv(args.cazymes, sep="\t", dtype=str)
        caz = set(caz_df["protein_id"] if "protein_id" in caz_df else caz_df["Gene ID"])  # pipeline / dbCAN overview
        n_caz = fam_df["protein"].isin(caz).sum()
        fam_df = fam_df[~fam_df["protein"].isin(caz)]
    label = fam_df.sort_values(["family", "length"], ascending=False).groupby("family", sort=False)["product"].first()
    bar = (fam_df.groupby("family")["plasmid"].nunique().rename("plasmids").reset_index()
           .sort_values("plasmids", ascending=False, kind="stable"))
    bar["product"] = bar["family"].map(label)
    bar = bar[~bar["product"].str.contains("conjuga", case=False, na=False)]
    bar.to_csv("signal_peptide_families.tsv", sep="\t", index=False)

    with open("signal_peptide_summary.tsv", "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"signal peptide proteins\t{len(df)}\n")
        out.write(f"signal peptide proteins in a Panaroo family\t{df['family'].notna().sum()}\n")
        out.write(f"families with a signal peptide protein\t{df['family'].nunique()}\n")
        out.write(f"plasmids with a signal peptide family (non-conjugative)\t{len(covered)}\n")
        out.write(f"clusters with a signal peptide family\t{n_cl}\n")
        out.write(f"clusters total\t{n_cl_total}\n")
        out.write(f"percent clusters\t{100 * n_cl / n_cl_total:.1f}\n")
        out.write(f"CAZymes removed for the family ranking\t{n_caz}\n")

    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = bar.head(args.top)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.bar(range(len(top)), top["plasmids"], width=0.6, color="#56B4E9")
    for i, v in enumerate(top["plasmids"]):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=6)
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top["product"], rotation=-45, ha="left", rotation_mode="anchor", fontsize=7)
    ax.set_xlabel(f"Protein family (top {args.top})")
    ax.set_ylabel("Number of accessory plasmids")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig("signal_peptide_families.png", dpi=300, bbox_inches="tight")
    fig.savefig("signal_peptide_families.svg", bbox_inches="tight")

    print(f"{len(df)} signal peptide proteins; {len(covered)} plasmids, {n_cl}/{n_cl_total} clusters "
          f"({100 * n_cl / n_cl_total:.1f}%)")


if __name__ == "__main__":
    main()
