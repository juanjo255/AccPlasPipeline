#!/usr/bin/env python3
"""
Per-plasmid summary and figures of the functional annotations of accessory plasmids:
toxin-antitoxin (TASmania), defense systems (DefenseFinder), biosynthetic gene
clusters (antiSMASH) and CAZymes (dbCAN, optional SignalP).

Clusters are the MobMess clusters of Step 4 (as in the publication notebooks).

Outputs:
  annotation_summary.tsv          one row per plasmid with counts and types (TA counts use
                                  genes with count = 1 in tasmania.ta.tsv)
  defense_systems.tsv, bgc_regions.tsv, cazymes_signalp.tsv  merged long tables
  Figures (PNG + SVG): ta_heatmap_per_plasmid, ta_heatmap_per_cluster, ta_operon_structures,
  defense_heatmap_per_cluster, bgc_types, cazyme_families
"""
import argparse
import csv
import glob
import os
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from cmcrameri import cm


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--plasmids", required=True, help="plasmid_summary.tsv (Step 9)")
    p.add_argument("--ta", required=True, help="tasmania.ta.tsv")
    p.add_argument("--defense", nargs="*", default=[], help="<plasmid>_defense_finder_systems.tsv files")
    p.add_argument("--antismash", required=True, help="antiSMASH output directory")
    p.add_argument("--cazymes", required=True, help="cazymes.tsv")
    p.add_argument("--signalp", help="optional SignalP 6 prediction_results.txt for the CAZymes")
    p.add_argument("--top", type=int, default=20, help="number of categories shown in heatmaps")
    return p.parse_args()


def save(fig, name):
    fig.savefig(f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def heatmap(data, name, title, xlabel, ylabel, cmap, cbar_label, row_labels=False, col_cluster=False):
    if data.empty or data.shape[0] < 2:
        print(f"{name}: not enough data, skipped", file=sys.stderr)
        return
    sns.set_theme(style="white")
    g = sns.clustermap(data, cmap=cmap, figsize=(9, 10), col_cluster=col_cluster and data.shape[1] > 1,
                       row_cluster=True, dendrogram_ratio=(0.10, 0.08 if col_cluster else 0.001),
                       cbar_pos=(1.02, 0.35, 0.02, 0.3), cbar_kws={"label": cbar_label},
                       xticklabels=True, yticklabels=row_labels)
    g.ax_heatmap.set_xlabel(xlabel)
    g.ax_heatmap.set_ylabel(ylabel)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    # Title above the dendrograms instead of on the heatmap axis
    g.figure.suptitle(title, y=1.02)
    g.savefig(f"{name}.png", dpi=300, bbox_inches="tight")
    g.savefig(f"{name}.svg", bbox_inches="tight")
    plt.close(g.figure)


def read_defense(paths):
    frames = []
    for path in paths:
        plasmid = os.path.basename(path).rsplit("_defense_finder_systems.tsv", 1)[0]
        df = pd.read_csv(path, sep="\t", dtype=str)
        if not df.empty:
            df.insert(0, "plasmid", plasmid)
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["plasmid", "type", "subtype"])


def read_antismash(outdir):
    # One row per (region, product); a region can have several products (hybrid BGCs)
    rows = []
    for path in sorted(glob.glob(os.path.join(outdir, "*.region*.gbk"))):
        m = re.match(r"(.+)\.region(\d+)\.gbk$", os.path.basename(path))
        record, region = m.group(1), int(m.group(2))
        in_region, products = False, []
        with open(path) as fh:
            for line in fh:
                if line.startswith("     ") and not line.startswith("      "):
                    in_region = line.split()[0] == "region"
                elif in_region:
                    pm = re.match(r'\s+/product="(.+)"', line)
                    if pm:
                        products.append(pm.group(1))
        for product in products or ["unknown"]:
            rows.append({"plasmid": record, "region": region, "product": product})
    return pd.DataFrame(rows, columns=["plasmid", "region", "product"])


def read_signalp(path):
    if not path:
        return {}
    pred = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            pred[f[0].split()[0]] = f[1]
    return pred


def main():
    args = parse_args()
    pls = pd.read_csv(args.plasmids, sep="\t", dtype=str).set_index("sequence")
    pls["length"] = pls["length"].astype(int)

    ta_all = pd.read_csv(args.ta, sep="\t", dtype=str, keep_default_na=False)
    # Counted TA genes: E < --evalue and not RepA/B/C (see tasmania_search.py)
    ta = ta_all[ta_all["count"] == "1"].copy()
    defense = read_defense(args.defense)
    bgc = read_antismash(args.antismash)
    caz = pd.read_csv(args.cazymes, sep="\t", dtype=str)
    sp = read_signalp(args.signalp)
    caz["signalp"] = caz["protein_id"].map(sp).fillna("not run" if not sp else "-")

    defense.to_csv("defense_systems.tsv", sep="\t", index=False)
    bgc.to_csv("bgc_regions.tsv", sep="\t", index=False)
    caz.to_csv("cazymes_signalp.tsv", sep="\t", index=False)

    def joined(df, col):
        return df.groupby("plasmid")[col].apply(lambda s: ";".join(sorted(set(s)))) if not df.empty else pd.Series(dtype=str)

    summ = pd.DataFrame(index=pls.index)
    summ["dataset"], summ["cluster"], summ["length"] = pls["dataset"], pls["cluster"], pls["length"]
    summ["n_ta_proteins"] = ta.groupby("plasmid").size()
    # Pseudo-operons holding both a toxin and an antitoxin (complete TA pairs)
    complete = ta[ta["operon_structure"].str.contains("T") & ta["operon_structure"].str.contains("A")]
    summ["n_ta_complete_operons"] = complete.drop_duplicates("operon_id").groupby("plasmid").size()
    summ["ta_operon_structures"] = joined(ta.drop_duplicates("operon_id"), "operon_structure")
    summ["n_defense_systems"] = defense.groupby("plasmid").size() if not defense.empty else 0
    summ["defense_types"] = joined(defense, "type")
    summ["n_bgc_regions"] = bgc.drop_duplicates(["plasmid", "region"]).groupby("plasmid").size() if not bgc.empty else 0
    summ["bgc_types"] = joined(bgc, "product")
    summ["n_cazymes"] = caz.groupby("plasmid").size()
    summ["cazyme_families"] = joined(caz, "family")
    for c in ("n_ta_proteins", "n_ta_complete_operons", "n_defense_systems", "n_bgc_regions", "n_cazymes"):
        summ[c] = summ[c].fillna(0).astype(int)
    summ = summ.fillna("-")
    summ.reset_index(names="plasmid").to_csv("annotation_summary.tsv", sep="\t", index=False)

    # TA heatmaps: Bakta product counts per plasmid, and their mean per MobMess cluster
    if not ta.empty:
        # Bakta product of the overlapping gene ('-' when Bakta has no matching gene)
        ta["product"] = ta["bakta_product"].replace({"": "no Bakta gene", "-": "no Bakta gene"})
        counts = ta.groupby(["plasmid", "product"]).size().unstack(fill_value=0)
        top = counts.sum().sort_values(ascending=False).index[: args.top]
        heatmap(counts[top], "ta_heatmap_per_plasmid", f"TA proteins per plasmid (top {len(top)} products)",
                "Protein annotation", f"Plasmids (n = {counts.shape[0]})", cm.berlin, "TA count per plasmid")
        all_counts = counts.reindex(pls.index, fill_value=0)
        mean = all_counts.groupby(pls["cluster"]).mean()
        top = mean.sum().sort_values(ascending=False).index[: args.top]
        mean = mean[top]
        mean = mean.loc[(mean != 0).any(axis=1)]
        heatmap(mean, "ta_heatmap_per_cluster", f"Mean TA proteins per MobMess cluster (top {len(top)})",
                "Protein annotation", f"Clusters (n = {mean.shape[0]})", cm.batlow, "Mean TA count",
                col_cluster=True)

    # TA operon structures: number of pseudo-operons of each structure (top N)
    if not ta.empty:
        ops = ta.drop_duplicates("operon_id")["operon_structure"].value_counts().head(args.top)
        sns.set_theme(style="ticks")
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(ops)), 4))
        ax.bar(ops.index, ops.values, color="#55a868")
        ax.set_title(f"TA pseudo-operon structures (top {len(ops)}; T toxin, A antitoxin, x other gene)")
        ax.set_xlabel("Operon structure")
        ax.set_ylabel("Pseudo-operons")
        ax.tick_params(axis="x", rotation=45)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        save(fig, "ta_operon_structures")

    # Defense systems: mean count per plasmid in each MobMess cluster
    if not defense.empty:
        counts = defense.groupby(["plasmid", "type"]).size().unstack(fill_value=0).reindex(pls.index, fill_value=0)
        mean = counts.groupby(pls["cluster"]).mean()
        top = mean.sum().sort_values(ascending=False).index[: args.top]
        mean = mean[top]
        mean = mean.loc[(mean != 0).any(axis=1)]
        heatmap(mean, "defense_heatmap_per_cluster", f"Mean defense systems per MobMess cluster (top {len(top)} types)",
                "Defense system type", f"Clusters (n = {mean.shape[0]})", "YlGnBu", "Mean count",
                row_labels=mean.shape[0] <= 40, col_cluster=True)

    # BGC types (each product of each region counted once)
    if not bgc.empty:
        counts = bgc["product"].value_counts()
        sns.set_theme(style="ticks")
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(counts)), 4))
        ax.bar(counts.index, counts.values, color="#4c72b0")
        ax.set_title(f"BGCs in accessory plasmids (n = {bgc['plasmid'].nunique()} plasmids)")
        ax.set_xlabel("BGC type")
        ax.set_ylabel("Count")
        ax.tick_params(axis="x", rotation=45)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        save(fig, "bgc_types")

    # CAZyme families: number of plasmids, stacked by SignalP prediction when available
    if not caz.empty:
        table = caz.drop_duplicates(["plasmid", "family", "signalp"]).groupby(["family", "signalp"]).size().unstack(fill_value=0)
        table = table.loc[table.sum(axis=1).sort_values(ascending=False).index]
        sns.set_theme(style="ticks")
        fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(table)), 4))
        table.plot(kind="bar", stacked=True, ax=ax, width=0.8,
                   color=sns.color_palette("colorblind", table.shape[1]))
        ax.set_title(f"CAZymes in accessory plasmids (n = {caz['plasmid'].nunique()} plasmids)")
        ax.set_xlabel("CAZyme family")
        ax.set_ylabel("Number of plasmids")
        ax.tick_params(axis="x", rotation=45)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
        ax.legend(title="SignalP", frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        save(fig, "cazyme_families")

    print(f"{len(summ)} plasmids: {int((summ['n_ta_proteins'] > 0).sum())} with TA, "
          f"{int((summ['n_defense_systems'] > 0).sum())} with defense systems, "
          f"{int((summ['n_bgc_regions'] > 0).sum())} with BGCs, {int((summ['n_cazymes'] > 0).sum())} with CAZymes")


if __name__ == "__main__":
    main()
