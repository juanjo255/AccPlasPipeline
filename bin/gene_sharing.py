#!/usr/bin/env python3
"""
Gene families shared between accessory plasmids and the main replicons (Panaroo).

Port of the publication's build_panaroo_gene_redundancy_tables.py. Panaroo is run on
the whole-genome annotations (one GFF per genome) plus one standalone GFF per
accessory plasmid, so every accessory-plasmid gene appears twice. To avoid double
counting, accessory-plasmid membership comes ONLY from the standalone plasmid columns
(column name = plasmid ID), and in the genome columns genes on accessory-plasmid
contigs are skipped; the other genome genes are assigned to Chromosome / pSymA /
pSymB through their contig tag. Refound genes are resolved through both
annotation_id and clustering_id of gene_data.csv ('_pseudo' suffix removed).

Outputs:
  family_group_membership.csv  Gene, Annotation, n_contigs_Chromosome, n_contigs_pSymA,
                               n_contigs_pSymB, n_contigs_AccPlasmid (all-zero families dropped)
  gene_sharing_summary.tsv     accessory-plasmid families shared with each main replicon
  gene_sharing_upset.png/.svg  UpSet plot of families per replicon group
  shared_families_cog.tsv / shared_families_cog.png/.svg  (with --eggnog) COG categories of
                               accessory-plasmid genes in families shared with each replicon
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)
GROUPS = ("Chromosome", "pSymA", "pSymB", "AccPlasmid")
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
    p.add_argument("--gpa", required=True, help="Panaroo gene_presence_absence.csv")
    p.add_argument("--gene-data", required=True, help="Panaroo gene_data.csv")
    p.add_argument("--tags", nargs="+", required=True,
                   help="contig tag tables: <sample>.replicon_tags.tsv (with header) or 'contig<TAB>tag' (no header)")
    p.add_argument("--plasmids", required=True, help="accessory plasmid IDs (one per line)")
    p.add_argument("--eggnog", help="eggNOG-mapper annotations of the accessory plasmid proteins (<locus>__<plasmid>)")
    p.add_argument("--gff-samples", help="optional 'gff_name<TAB>sample': genome GFF contigs without the _<sample> "
                                         "suffix are matched to the tag table after adding it")
    return p.parse_args()


def norm(token):
    return token[:-7] if token.endswith("_pseudo") else token


def main():
    args = parse_args()
    contig2replicon = {}
    for path in args.tags:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or f[0] == "contig":
                    continue
                rep = f[1].removesuffix("_fragment")  # secondary pieces count for their replicon
                if rep in ("Chromosome", "pSymA", "pSymB"):
                    contig2replicon[f[0]] = rep
    with open(args.plasmids) as fh:
        acc_contigs = {l.strip() for l in fh if l.strip()}

    sample_of_gff = {}
    if args.gff_samples:
        with open(args.gff_samples) as fh:
            sample_of_gff = dict(l.rstrip("\n").split("\t")[:2] for l in fh if l.strip())

    id2scaffold = {}
    with open(args.gene_data) as fh:
        r = csv.reader(fh)
        h = next(r)
        i_s, i_c, i_a, i_g = h.index("scaffold_name"), h.index("clustering_id"), h.index("annotation_id"), h.index("gff_file")
        for row in r:
            scaffold = row[i_s]
            smp = sample_of_gff.get(row[i_g])
            if smp and scaffold not in contig2replicon and scaffold not in acc_contigs \
                    and not scaffold.endswith("_" + smp):
                scaffold = f"{scaffold}_{smp}"  # same renaming as the contig length filter
            id2scaffold[row[i_c]] = scaffold
            if row[i_a]:
                id2scaffold[row[i_a]] = scaffold

    family_contigs = defaultdict(lambda: defaultdict(set))
    family_annot, acc_gene_family = {}, []
    n_unresolved = 0
    with open(args.gpa) as fh:
        r = csv.reader(fh)
        header = next(r)
        sample_cols = header[3:]
        is_acc = [c in acc_contigs for c in sample_cols]
        if not any(is_acc):
            sys.exit("ERROR: no Panaroo sample column matches an accessory plasmid ID (standalone plasmid GFFs missing?)")
        for row in r:
            gene, annotation = row[0], row[2]
            family_annot[gene] = annotation
            groups = family_contigs[gene]
            for i, cell in enumerate(row[3:]):
                if not cell:
                    continue
                if is_acc[i]:
                    groups["AccPlasmid"].add(sample_cols[i])
                    for token in cell.split(";"):
                        acc_gene_family.append((sample_cols[i], token.strip(), gene))
                    continue
                for token in cell.split(";"):
                    scaffold = id2scaffold.get(norm(token))
                    if scaffold is None:
                        n_unresolved += 1
                        continue
                    if scaffold in acc_contigs:
                        continue
                    rep = contig2replicon.get(scaffold)
                    if rep is not None:
                        groups[rep].add(scaffold)

    sets = {g: set() for g in GROUPS}
    with open("family_group_membership.csv", "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["Gene", "Annotation"] + [f"n_contigs_{g}" for g in GROUPS])
        for gene, annotation in family_annot.items():
            counts = [len(family_contigs.get(gene, {}).get(g, ())) for g in GROUPS]
            if sum(counts) == 0:
                continue
            w.writerow([gene, annotation] + counts)
            for g, n in zip(GROUPS, counts):
                if n:
                    sets[g].add(gene)

    pls = sets["AccPlasmid"]
    any_main = pls & (sets["Chromosome"] | sets["pSymA"] | sets["pSymB"])
    with open("gene_sharing_summary.tsv", "w") as out:
        out.write("comparison\tfamilies\tpercent_of_plasmid_families\n")
        out.write(f"accessory plasmid families\t{len(pls)}\t100.0\n")
        for g in ("Chromosome", "pSymA", "pSymB"):
            n = len(pls & sets[g])
            out.write(f"shared with {g}\t{n}\t{100 * n / len(pls):.1f}\n")
        out.write(f"shared with any main replicon\t{len(any_main)}\t{100 * len(any_main) / len(pls):.1f}\n")
        out.write(f"accessory plasmids only\t{len(pls) - len(any_main)}\t{100 * (len(pls) - len(any_main)) / len(pls):.1f}\n")

    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import upsetplot

    contents = {"Chromosome": sets["Chromosome"], "pSymA": sets["pSymA"], "pSymB": sets["pSymB"],
                "Acc. Plasmids": sets["AccPlasmid"]}
    up = upsetplot.UpSet(upsetplot.from_contents(contents), sort_by="cardinality", sort_categories_by="cardinality",
                         show_counts="%d", facecolor="#56B4E9", element_size=None)
    up.style_subsets(present="Acc. Plasmids", facecolor="#E69F00", label="Involves Acc. Plasmids")
    fig = plt.figure(figsize=(12, 5))
    axes = up.plot(fig)
    axes["intersections"].set_ylabel("Gene families")
    fig.savefig("gene_sharing_upset.png", dpi=300, bbox_inches="tight")
    fig.savefig("gene_sharing_upset.svg", bbox_inches="tight")
    plt.close(fig)

    if args.eggnog:
        # COG letters of accessory-plasmid genes in families shared with each replicon
        with open(args.eggnog) as fh:
            rows = [l for l in fh if not l.startswith("##")]
        cog_of = {}
        for r in csv.DictReader(rows, delimiter="\t"):
            cog_of[r["#query"].split("__", 1)[0]] = r["COG_category"] or "-"
        recs = []
        for contig, locus, family in acc_gene_family:
            if locus in cog_of:
                for letter in cog_of[locus]:
                    if letter not in ("-", "S"):
                        for g in ("Chromosome", "pSymA", "pSymB"):
                            if family in sets[g]:
                                recs.append((f"shared with {g}", contig, family, letter, COG.get(letter, letter)))
        df = pd.DataFrame(recs, columns=["comparison", "plasmid", "family", "cog", "description"]).drop_duplicates()
        # one count per plasmid per family per COG letter, as in the publication notebook
        tab = (df.groupby(["comparison", "family", "description"]).size().reset_index(name="plasmids")
               .groupby(["description", "comparison"])["plasmids"].max().unstack(fill_value=0))
        tab.to_csv("shared_families_cog.tsv", sep="\t")
        if not tab.empty:
            tab = tab.loc[tab.max(axis=1).sort_values(ascending=False).index]
            fig, ax = plt.subplots(figsize=(10, 4.5))
            tab.plot(kind="bar", ax=ax, width=0.8)
            ax.set_ylabel("Max accessory plasmids per gene family")
            ax.set_xlabel("COG category")
            ax.tick_params(axis="x", rotation=60, labelsize=7)
            for lab in ax.get_xticklabels():
                lab.set_ha("right")
            ax.legend(frameon=False)
            ax.spines[["top", "right"]].set_visible(False)
            fig.savefig("shared_families_cog.png", dpi=300, bbox_inches="tight")
            fig.savefig("shared_families_cog.svg", bbox_inches="tight")
            plt.close(fig)

    print(f"{len(family_annot)} families; accessory plasmid families {len(pls)}, shared with a main replicon "
          f"{len(any_main)} ({100 * len(any_main) / len(pls):.1f}%); unresolved genome-side ids {n_unresolved}")


if __name__ == "__main__":
    main()
