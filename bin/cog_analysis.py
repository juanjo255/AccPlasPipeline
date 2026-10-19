#!/usr/bin/env python3
"""
COG content of accessory plasmids from eggNOG-mapper (publication notebook methods).

  - COG category counts: every letter of a multi-letter category (e.g. 'KL')
    counts once; '-' is 'No category assigned'.
  - Gene content: presence/absence of the root eggNOG orthologous group (first
    entry of eggNOG_OGs, e.g. 'COG1028@1|root' -> 'COG1028') per plasmid;
    Jaccard distances, PCoA (classical scaling) and metric MDS (random_state 0).
  - Enrichment: for each MobMess cluster with more than --min-cluster members, a
    one-sided Fisher exact test (greater) per COG letter of the number of plasmids
    carrying it, cluster vs all other plasmids; log2FC of the proportions with a
    1e-10 pseudocount; Benjamini-Hochberg over all tests. Significant: FDR < 0.05
    and log2FC > 0.2.
  - Cluster composition heatmap: COG letter counts per cluster (clusters with at
    least --min-cluster members; 'S' and '-' left out), row z-scores.

Protein IDs are <locus>__<plasmid> (Step 10 merged annotation).

Outputs: cog_counts.tsv, og_presence_absence.tsv, cog_ordination.tsv,
cog_enrichment.tsv (all tests), cog_enrichment_significant.tsv, figures (PNG + SVG):
cog_counts, cog_pcoa, cog_mds, cog_enrichment_dotplot, cog_cluster_heatmap
"""
import argparse
import os
import sys
import warnings

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact
from sklearn.manifold import MDS
from statsmodels.stats.multitest import multipletests

COG = {
    "J": "Translation, ribosomal structure and biogenesis", "A": "RNA processing and modification",
    "K": "Transcription", "L": "Replication, recombination and repair", "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning", "Y": "Nuclear structure",
    "V": "Defense mechanisms", "T": "Signal transduction mechanisms", "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility", "Z": "Cytoskeleton", "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones", "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism", "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism", "R": "General function prediction only",
    "S": "Function unknown", "-": "No category assigned",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotations", required=True, help="eggNOG-mapper *.emapper.annotations")
    p.add_argument("--plasmids", required=True, help="plasmid_summary.tsv (sequence, cluster)")
    p.add_argument("--min-cluster", type=int, default=3, help="cluster size limit (enrichment: > N, heatmap: >= N)")
    p.add_argument("--top-clusters", type=int, default=10, help="clusters coloured in the ordination plots")
    return p.parse_args()


def read_emapper(path):
    # Comment lines start with '##'; the header line starts with '#query'
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("##")]
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str, keep_default_na=False)
    df = df.rename(columns={"#query": "query"})
    df["plasmid"] = df["query"].str.split("__", n=1).str[1]
    df["locus_tag"] = df["query"].str.split("__", n=1).str[0]
    if df["plasmid"].isna().any():
        sys.exit("ERROR: eggNOG query IDs must be <locus>__<plasmid>")
    return df


def save(fig, name):
    fig.savefig(f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def pcoa(dist, n=2):
    # Classical multidimensional scaling (double-centred eigendecomposition)
    m = dist.shape[0]
    j = np.eye(m) - np.ones((m, m)) / m
    b = -0.5 * j @ (dist ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals > 0
    coords = vecs[:, :n] * np.sqrt(np.clip(vals[:n], 0, None))
    explained = vals[:n] / vals[pos].sum() * 100
    return coords, explained


def main():
    args = parse_args()
    ann = read_emapper(args.annotations)
    pls = pd.read_csv(args.plasmids, sep="\t", dtype=str).set_index("sequence")
    cluster_of = pls["cluster"]

    # 1. COG category counts (each letter once per protein)
    letters = ann["COG_category"].replace("", "-").apply(list).explode()
    counts = letters.value_counts().rename_axis("cog").reset_index(name="count")
    counts["description"] = counts["cog"].map(COG)
    counts.to_csv("cog_counts.tsv", sep="\t", index=False)
    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(counts["description"], counts["count"], color="#4c72b0")
    ax.set_title(f"COG categories of accessory plasmid proteins (n = {ann['plasmid'].nunique()} plasmids)")
    ax.set_ylabel("Proteins")
    ax.tick_params(axis="x", rotation=60, labelsize=7)
    for lab in ax.get_xticklabels():
        lab.set_ha("right")
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "cog_counts")

    # COG letter counts per plasmid (enrichment uses presence of each letter)
    per_pls = (ann.assign(cog=ann["COG_category"].replace("", "-").apply(list)).explode("cog")
               .groupby(["plasmid", "cog"]).size().unstack(fill_value=0))

    # 2. Root orthologous groups -> Jaccard -> PCoA / MDS
    og = ann["eggNOG_OGs"].where(~ann["eggNOG_OGs"].isin(["", "-"]))
    ann["OG"] = og.str.split(",").str[0].str.split("@").str[0]
    og_pa = (ann.dropna(subset=["OG"]).groupby(["plasmid", "OG"]).size().unstack(fill_value=0) > 0).astype(int)
    og_pa.to_csv("og_presence_absence.tsv", sep="\t")
    dist = squareform(pdist(og_pa.to_numpy().astype(bool), metric="jaccard"))
    pc, expl = pcoa(dist)
    mds = MDS(n_components=2, random_state=0, dissimilarity="precomputed").fit_transform(dist)
    ordn = pd.DataFrame({"plasmid": og_pa.index, "cluster": cluster_of.reindex(og_pa.index).values,
                         "PCoA_1": pc[:, 0], "PCoA_2": pc[:, 1], "MDS_1": mds[:, 0], "MDS_2": mds[:, 1]})
    ordn.to_csv("cog_ordination.tsv", sep="\t", index=False)

    sizes = ordn["cluster"].value_counts()
    keep = sizes[sizes >= 5].index
    plot_df = ordn[ordn["cluster"].isin(keep)].copy()
    top = plot_df["cluster"].value_counts().index[: args.top_clusters]
    plot_df["colour"] = plot_df["cluster"].where(plot_df["cluster"].isin(top), "Other")
    pal = dict(zip(top, sns.color_palette("tab10", len(top))))
    pal["Other"] = "#D7D7D7"
    for x, y, name, xl, yl in (("PCoA_1", "PCoA_2", "cog_pcoa", f"PCo1 ({expl[0]:.1f}%)", f"PCo2 ({expl[1]:.1f}%)"),
                               ("MDS_1", "MDS_2", "cog_mds", "MDS1", "MDS2")):
        fig, ax = plt.subplots(figsize=(6, 5))
        for c in ["Other"] + list(top):
            sub = plot_df[plot_df["colour"] == c]
            ax.scatter(sub[x], sub[y], s=10, color=pal[c], label=f"cluster {c}" if c != "Other" else "Other",
                       linewidths=0, alpha=0.9)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(f"Gene content (Jaccard on eggNOG OGs), top {len(top)} MobMess clusters (n >= 5)", fontsize=9)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=2, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        save(fig, name)

    # 3. COG enrichment per MobMess cluster vs the rest (notebook compute_COG_enrichment)
    universe = list(og_pa.index)  # plasmids in the gene-content analysis
    cl_sizes = cluster_of.value_counts()
    valid = set(cl_sizes[cl_sizes > args.min_cluster].index)
    records = []
    in_clusters = ordn[ordn["cluster"].isin(valid)]
    for cid in in_clusters["cluster"].unique():
        in_idx = list(in_clusters.loc[in_clusters["cluster"] == cid, "plasmid"])
        out_idx = [p for p in universe if p not in set(in_idx)]
        n_in, n_out = len(in_idx), len(out_idx)
        for cog in per_pls.columns:
            a = int((per_pls.reindex(in_idx)[cog].fillna(0) > 0).sum())
            c = int((per_pls.reindex(out_idx)[cog].fillna(0) > 0).sum())
            b, d = n_in - a, n_out - c
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            pin, pout = (a / n_in if n_in else 0), (c / n_out if n_out else 0)
            records.append({"cluster": cid, "cog": cog, "description": COG.get(cog, cog), "p_value": p,
                            "log2fc": np.log2((pin + 1e-10) / (pout + 1e-10)), "a": a, "b": b, "c": c, "d": d})
    enr = pd.DataFrame(records)
    if not enr.empty:
        enr["p_adj"] = multipletests(enr["p_value"], method="fdr_bh")[1]
    enr.to_csv("cog_enrichment.tsv", sep="\t", index=False)
    sig = enr[(enr["p_adj"] < 0.05) & (enr["log2fc"] > 0.2)].sort_values(["cluster", "p_adj"]) if not enr.empty else enr
    sig.to_csv("cog_enrichment_significant.tsv", sep="\t", index=False)

    dot = sig[~sig["cog"].isin(["S", "-"])] if not sig.empty else sig
    if not dot.empty:
        dot = dot.assign(mlog=-np.log10(dot["p_adj"].clip(lower=1e-30)))
        clusters = sorted(dot["cluster"].unique(), key=lambda s: int(s) if str(s).isdigit() else s)
        cats = sorted(dot["description"].unique())
        fig, ax = plt.subplots(figsize=(max(4, 0.45 * len(clusters) + 3), max(3, 0.35 * len(cats) + 1)))
        sc = ax.scatter([clusters.index(c) for c in dot["cluster"]], [cats.index(c) for c in dot["description"]],
                        s=20 + 25 * dot["mlog"], c=dot["log2fc"], cmap="Reds", edgecolors="grey", linewidths=0.3)
        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters, rotation=90)
        ax.set_yticks(range(len(cats)))
        ax.set_yticklabels(cats)
        ax.set_xlabel("MobMess cluster")
        fig.colorbar(sc, ax=ax, label="log2FC")
        for v in (1, 5, 10):
            ax.scatter([], [], s=20 + 25 * v, color="lightgrey", edgecolors="grey", label=f"-log10 FDR = {v}")
        ax.legend(frameon=False, bbox_to_anchor=(1.25, 1), loc="upper left", fontsize=8)
        ax.set_title("COG categories enriched per cluster (FDR < 0.05, log2FC > 0.2)", fontsize=9)
        save(fig, "cog_enrichment_dotplot")

    # 4. COG composition per cluster (clusters with >= min_cluster members), row z-score
    cl_ok = set(cl_sizes[cl_sizes >= args.min_cluster].index)
    long = ann.assign(cluster=ann["plasmid"].map(cluster_of), cog=ann["COG_category"].replace("", "-").apply(list))
    long = long[long["cluster"].isin(cl_ok)].explode("cog")
    comp = pd.crosstab(long["cluster"], long["cog"]).drop(columns=["S", "-"], errors="ignore").rename(columns=COG)
    comp = comp.loc[:, comp.sum() > 0]
    if comp.shape[0] > 1 and comp.shape[1] > 1:
        g = sns.clustermap(comp, cmap="viridis", z_score=0, figsize=(10, 12), cbar_pos=(1.02, 0.8, 0.02, 0.15))
        g.ax_cbar.set_ylabel("row z-score")
        g.ax_heatmap.set_xlabel("COG category")
        g.ax_heatmap.set_ylabel(f"MobMess clusters (n = {comp.shape[0]})")
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=6)
        g.savefig("cog_cluster_heatmap.png", dpi=300, bbox_inches="tight")
        g.savefig("cog_cluster_heatmap.svg", bbox_inches="tight")
        plt.close(g.figure)

    print(f"{ann['plasmid'].nunique()} plasmids, {len(ann)} annotated proteins; {og_pa.shape[1]} root OGs; "
          f"{len(valid)} clusters tested, {len(enr)} tests, {len(sig)} significant")


if __name__ == "__main__":
    main()
