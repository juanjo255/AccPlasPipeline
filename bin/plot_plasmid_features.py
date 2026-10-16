#!/usr/bin/env python3
"""
Summary table and figures of the accessory plasmids (kept in the plasmid-marker step).

Figures (PNG and SVG), adapted from the publication notebooks:
  features_stackbar           % of plasmids per relaxase type, MPF type, predicted
                              mobility, host range (MOB-typer) and rep genes (Bakta)
  rep_proteins_heatmap        replication/partition proteins per plasmid, rows
                              clustered (Hamming, average) and coloured by dataset
  replicon_groups_by_dataset  plasmids per replicon group, stacked by dataset
  replicon_multiplicity_length  length of plasmids with single vs duplicated
                              repA/B/C (Mann-Whitney U test)
  replicon_cooccurrence_upset UpSet plot of replicon groups carried together by
                              the same strain (sample plasmids with a replicon group)
  plasmids_per_strain         strains by number of accessory plasmids
                              (strains with at least one)

Plasmids are coloured by the samplesheet `dataset` column; PLSDB plasmids are 'Database'.

Table: plasmid_summary.tsv, one row per accessory plasmid.
"""
import argparse
import csv
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
import upsetplot

DATABASE = "Database"
MISSING_LABEL, OTHER_LABEL = "Missing", "Other"
MISSING_COLOR, OTHER_COLOR, DATABASE_COLOR = "#D2D2D2", "#0ccf8b", "#aaaaaa"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--features", required=True, help="plasmid_features.features.tsv")
    p.add_argument("--rep-proteins", required=True, help="plasmid_features.rep_proteins.tsv")
    p.add_argument("--mobtyper", required=True, help="MOB-typer <prefix>_mob_type.tsv")
    p.add_argument("--repabc", required=True, help="repabc.summary.tsv")
    p.add_argument("--replicon-groups", required=True, help="repA.replicon_groups.tsv")
    p.add_argument("--datasets", required=True, help="contig <tab> sample <tab> dataset")
    return p.parse_args()


def save(fig, name):
    fig.savefig(f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def dataset_palette(datasets):
    # Colour-blind palette per sample dataset; PLSDB always grey and last
    names = sorted(d for d in set(datasets) if d != DATABASE)
    pal = dict(zip(names, sns.color_palette("colorblind", max(len(names), 1))))
    pal[DATABASE] = DATABASE_COLOR
    return pal


def build_table(args):
    feat = pd.read_csv(args.features, sep="\t", dtype=str)
    feat = feat[feat["keep"] == "1"].copy()
    feat["length"] = feat["length"].astype(int)

    ds = pd.read_csv(args.datasets, sep="\t", dtype=str).drop_duplicates("contig").set_index("contig")
    feat["sample"] = feat["sequence"].map(ds["sample"])
    feat["dataset"] = feat["sequence"].map(ds["dataset"])
    feat.loc[feat["source"] == "plsdb", ["sample", "dataset"]] = ["-", DATABASE]
    missing = feat["dataset"].isna()
    if missing.any():
        sys.exit(f"ERROR: no dataset for: {', '.join(feat.loc[missing, 'sequence'][:5])}")

    mob = pd.read_csv(args.mobtyper, sep="\t", dtype=str).set_index("sample_id")
    feat["host_range"] = feat["sequence"].map(mob["observed_host_range_ncbi_name"]) if "observed_host_range_ncbi_name" in mob else "-"

    rep = pd.read_csv(args.repabc, sep="\t", dtype=str).set_index("sequence")
    feat["rep_genes"] = feat["sequence"].map(rep["rep_genes"]).fillna("No RepABC")
    feat["complete_repabc"] = feat["sequence"].map(rep["complete_operon"]).fillna("0")

    groups = pd.read_csv(args.replicon_groups, sep="\t", dtype=str).set_index("sequence")
    feat["replicon_group"] = feat["sequence"].map(groups["replicon_group"]).fillna("-")

    cols = ["sequence", "sample", "dataset", "source", "length", "circular", "cluster", "cluster_type",
            "replicon_group", "rep_genes", "complete_repabc", "rep_partition_proteins",
            "mobtyper_relaxase", "mobtyper_mpf", "mobtyper_mobility", "host_range",
            "genomad_plasmid_score", "genomad_mob", "rel_coverage", "review_flag", "reason"]
    feat[cols].sort_values(["dataset", "sequence"]).to_csv("plasmid_summary.tsv", sep="\t", index=False)
    return feat


def features_stackbar(df):
    panels = [("mobtyper_relaxase", "Relaxase type", None), ("mobtyper_mpf", "MPF type", None),
              ("mobtyper_mobility", "Predicted mobility", None), ("host_range", "Host range", 8),
              ("rep_genes", "Rep genes", None)]
    sns.set_theme(style="ticks", font_scale=1.1)
    fig, axes = plt.subplots(len(panels), 1, figsize=(15, 7.5), sharex=True)
    for ax, (col, title, top_n) in zip(axes, panels):
        values = df[col].replace({"-": MISSING_LABEL, "NA": MISSING_LABEL}).fillna(MISSING_LABEL)
        if top_n is not None:
            keep = values.value_counts().index[:top_n]
            values = values.where(values.isin(keep), OTHER_LABEL)
        pct = values.value_counts(normalize=True).mul(100).to_frame().T
        cats = pct.columns.tolist()
        regular = [c for c in cats if c not in (MISSING_LABEL, OTHER_LABEL)]
        pal = dict(zip(regular, sns.color_palette("colorblind", max(len(regular), 1))))
        pal.update({MISSING_LABEL: MISSING_COLOR, OTHER_LABEL: OTHER_COLOR})
        pct.plot(kind="barh", stacked=True, ax=ax, width=0.5, color=[pal[c] for c in cats])
        ax.set_yticks([])
        ax.set_xlim(0, 100)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10, ncol=2, title=title)
    axes[-1].set_xlabel(f"Percentage of plasmids (%, n = {len(df)})")
    fig.tight_layout()
    save(fig, "features_stackbar")


def rep_proteins_heatmap(df, rep_path, palette):
    rp = pd.read_csv(rep_path, sep="\t", dtype=str)
    rp = rp[rp["Sequence Id"].isin(df["sequence"])].copy()
    if rp.empty:
        print("No replication/partition proteins: heatmap skipped", file=sys.stderr)
        return
    rp["protein"] = rp["Gene"].fillna("").replace("", OTHER_LABEL)
    counts = rp.groupby(["Sequence Id", "protein"]).size().unstack(fill_value=0).astype(int)
    datasets = df.set_index("sequence")["dataset"].reindex(counts.index)
    row_colors = datasets.map(palette).rename("Dataset")

    max_val = int(counts.values.max())
    ramp = ["#c3c5c5", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#184f95"]
    cmap = ListedColormap((ramp * (max_val // len(ramp) + 1))[: max_val + 1])
    bounds = list(range(0, max_val + 2))
    norm = BoundaryNorm(bounds, cmap.N)

    sns.set_theme(style="white", font_scale=1.05)
    g = sns.clustermap(counts, cmap=cmap, norm=norm, col_cluster=False, row_cluster=len(counts) > 1,
                       method="average", metric="hamming", row_colors=row_colors, yticklabels=False,
                       xticklabels=True, linewidths=0, figsize=(6, 9), dendrogram_ratio=(0.12, 0.06),
                       cbar_pos=(0.34, 0.98, 0.22, 0.02), cbar_kws={"orientation": "horizontal"})
    g.ax_heatmap.collections[0].set_edgecolor("face")
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel(f"Plasmids (n = {counts.shape[0]})")
    g.cax.set_xticks([b + 0.5 for b in bounds[:-1]])
    g.cax.set_xticklabels([str(b) for b in bounds[:-1]], fontsize=9)
    g.cax.set_title("Count", fontsize=10, pad=6)
    used = [d for d in palette if d in set(datasets)]
    g.ax_heatmap.legend(handles=[Patch(facecolor=palette[d], label=d) for d in used], title="Dataset",
                        bbox_to_anchor=(1.28, 1.0), loc="upper left", frameon=False, fontsize=9)
    g.figure.suptitle("Replication-related proteins per plasmid", y=1.05)
    g.savefig("rep_proteins_heatmap.png", dpi=300, bbox_inches="tight")
    g.savefig("rep_proteins_heatmap.svg", bbox_inches="tight")
    plt.close(g.figure)
    return counts


def replicon_groups_by_dataset(df, palette):
    typed = df[df["replicon_group"] != "-"]
    if typed.empty:
        print("No replicon groups: plot skipped", file=sys.stderr)
        return
    table = typed.groupby(["replicon_group", "dataset"]).size().unstack(fill_value=0)
    table = table.loc[table.sum(axis=1).sort_values(ascending=False).index]
    order = [d for d in palette if d in table.columns]
    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=(8, 4))
    table[order].plot(kind="bar", stacked=True, ax=ax, width=0.8, color=[palette[d] for d in order])
    ax.set_title("Replicon groups")
    ax.set_xlabel("Replicon group")
    ax.set_ylabel("Number of plasmids")
    ax.tick_params(axis="x", rotation=45)
    for lab in ax.get_xticklabels():
        lab.set_ha("right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="Dataset", frameon=False)
    fig.tight_layout()
    save(fig, "replicon_groups_by_dataset")


def replicon_multiplicity_length(df, counts):
    if counts is None:
        return
    repabc = counts[[c for c in ("repA", "repB", "repC") if c in counts.columns]]
    repabc = repabc[repabc.sum(axis=1) > 0]
    if repabc.empty:
        return
    cls = repabc.max(axis=1).map(lambda m: "Multi replicon" if m > 1 else "Single replicon")
    data = pd.DataFrame({"Replicon": cls, "Length": df.set_index("sequence")["length"].reindex(cls.index)})
    order = ["Single replicon", "Multi replicon"]
    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=(4, 4.5))
    sns.violinplot(data=data, x="Replicon", y="Length", order=order, cut=0, ax=ax, color="#9ec5f4")
    title = "Plasmid length vs RepABC multiplicity"
    groups = [data.loc[data["Replicon"] == o, "Length"] for o in order]
    if all(len(g) > 0 for g in groups):
        u, p = mannwhitneyu(*groups, alternative="two-sided")
        title += f"\nMann-Whitney U = {u:.0f}, p = {p:.3g}"
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{o}\n(n = {(data['Replicon'] == o).sum()})" for o in order])
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("Length (bp)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "replicon_multiplicity_length")


def replicon_cooccurrence_upset(df):
    typed = df[(df["source"] == "sample") & (df["replicon_group"] != "-")]
    if typed["sample"].nunique() < 2:
        print("Fewer than two strains with replicon groups: UpSet plot skipped", file=sys.stderr)
        return
    # Each strain is one element; its memberships are the replicon groups it carries
    memberships = typed.groupby("sample")["replicon_group"].apply(lambda g: sorted(set(g)))
    data = upsetplot.from_memberships(memberships.tolist())
    fig = plt.figure(figsize=(12, 6))
    axes = upsetplot.UpSet(data, subset_size="count", show_counts="{:,}", sort_by="cardinality",
                           facecolor="#4daf4a", element_size=None).plot(fig=fig)
    axes["intersections"].set_title("Replicon group co-occurrence across strains")
    axes["intersections"].set_ylabel("Strains")
    axes["totals"].set_xlabel("Strains")
    for ax in axes.values():
        if ax is not None:
            ax.grid(False)
    save(fig, "replicon_cooccurrence_upset")


def plasmids_per_strain(df):
    per_strain = df[df["source"] == "sample"].groupby("sample").size()
    if per_strain.empty:
        return
    counts = per_strain.value_counts().sort_index()
    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=(6, 3))
    bars = ax.barh([str(k) for k in counts.index], counts.values, color="#c44e52")
    ax.bar_label(bars, padding=3)
    ax.invert_yaxis()
    ax.set_title(f"Accessory plasmids per strain ({len(per_strain)} strains)")
    ax.set_xlabel("Number of strains")
    ax.set_ylabel("Number of plasmids")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "plasmids_per_strain")


def main():
    args = parse_args()
    df = build_table(args)
    palette = dataset_palette(df["dataset"])
    features_stackbar(df)
    counts = rep_proteins_heatmap(df, args.rep_proteins, palette)
    replicon_groups_by_dataset(df, palette)
    replicon_multiplicity_length(df, counts)
    replicon_cooccurrence_upset(df)
    plasmids_per_strain(df)
    print(f"{len(df)} plasmids plotted; datasets: "
          + ", ".join(f"{d} ({n})" for d, n in df["dataset"].value_counts().items()))


if __name__ == "__main__":
    main()
