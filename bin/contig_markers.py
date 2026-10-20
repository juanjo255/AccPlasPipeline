#!/usr/bin/env python3
"""
Plasmid markers of every filtered sample contig, and linear contigs with RepABC genes
that the MobMess cluster filter dropped (for user review).

One row per sample contig (contig_markers.tsv):
  sample, dataset, length, circular (Flye), assembly_contigs / assembly_circular_contigs
  (whole assembly), repA / repB / repC genes (RepABC HMMs), rep_genes and complete_operon
  (publication operon rules), geNomad plasmid score / hallmarks / conjugation genes,
  rel_coverage (depth / chromosome depth), MobMess cluster and cluster type, kept by the
  MobMess filter, replicon tag (Chromosome / pSymA / pSymB / Other, with a reference genome).

With --graph (Flye assembly-graph context, assembly_graph.py): the other contigs in the
contig's graph component (named by replicon tag when available), whether the component is
isolated, and the repeat edges on its path.

Review (linear_marker_review.tsv), linear contigs not tagged as a main replicon (or a
fragment of one, e.g. pSymA_fragment) that either
  A. carry >= 1 RepABC gene and were dropped by the MobMess filter, or
  B. carry a plasmid marker (>= 1 RepABC gene, or a geNomad plasmid call with >= 1 hallmark gene
     or, without hallmarks, a score >= --genomad-min-score) and their graph
     component holds no main replicon: alone, or only with other non-main contigs (with replicon
     tags; without tags: alone). Probably a complete replicon that Flye could not circularize,
     e.g. because of repeats. Kept or dropped. No length threshold: small plasmids are real too
     (e.g. the 7.5 kb contig_8_TG17).
Graph connections are reported, not used to exclude: a contig connected to the chromosome is
likely a fragment or an integrated element; a connection to pSymA/pSymB is ambiguous.
The default selection is not changed: these contigs are only reported, with the reasons.
"""
import argparse
import csv
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)
MAIN = {"Chromosome", "pSymA", "pSymB",
        "Chromosome_fragment", "pSymA_fragment", "pSymB_fragment"}  # fragments: secondary replicon pieces


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datasets", required=True, help="contig, sample, dataset (with header)")
    p.add_argument("--contig-info", nargs="+", required=True, help="<sample>.contig_info.tsv (contig, circular, length)")
    p.add_argument("--repabc", nargs="+", required=True, help="<sample>.repabc.tsv (RepABC HMM genes, Bakta layout)")
    p.add_argument("--operons", nargs="+", required=True, help="<sample>.summary.tsv (repabc_operons.py)")
    p.add_argument("--assembly-stats", nargs="+", required=True, help="<sample>.assembly_stats.tsv")
    p.add_argument("--genomad", nargs="*", default=[], help="geNomad <sample>_plasmid_summary.tsv")
    p.add_argument("--coverage", nargs="*", default=[], help="<sample>.coverage.tsv (FLYE_CIRCULARITY)")
    p.add_argument("--clusters", required=True, help="MobMess <prefix>_clusters.txt (all clusters)")
    p.add_argument("--kept", required=True, help="IDs kept by the MobMess cluster filter")
    p.add_argument("--tags", nargs="*", default=[], help="<sample>.replicon_tags.tsv (optional)")
    p.add_argument("--graph", nargs="*", default=[], help="<sample>.graph.tsv (assembly_graph.py, optional)")
    p.add_argument("--genomad-min-score", type=float, default=0.9,
                   help="criterion B: a geNomad plasmid call without hallmark genes is a marker from this score")
    return p.parse_args()


def read_tsv(paths):
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.DictReader((l for l in fh if not l.startswith("#") or l.startswith("#Sequence")),
                                      delimiter="\t"):
                yield row


def main():
    args = parse_args()
    meta = {r["contig"]: r for r in read_tsv([args.datasets])}
    info = {r["contig"]: r for r in read_tsv(args.contig_info)}
    stats = {r["sample"]: r for r in read_tsv(args.assembly_stats)}

    rep = defaultdict(lambda: {"repA": 0, "repB": 0, "repC": 0})
    for r in read_tsv(args.repabc):
        rep[r["#Sequence Id"]][r["Gene"]] += 1
    operon = {r["sequence"]: r for r in read_tsv(args.operons)}
    genomad = {r["seq_name"]: r for r in read_tsv(args.genomad)}
    cov = {r["contig"]: r for r in read_tsv(args.coverage)}
    tags = {}
    for r in read_tsv(args.tags):
        tags[r["contig"]] = r["tag"]

    graph = {r["contig"]: r for r in read_tsv(args.graph)}
    clusters = {r["contig"]: r for r in read_tsv([args.clusters])}
    members = defaultdict(list)
    for c, r in clusters.items():
        members[r["cluster"]].append(r)
    with open(args.kept) as fh:
        kept = {l.strip() for l in fh if l.strip()}

    cols = ["contig", "sample", "dataset", "length", "circular", "assembly_contigs", "assembly_circular_contigs",
            "repA", "repB", "repC", "rep_genes", "complete_operon", "genomad_plasmid_score", "genomad_hallmarks",
            "genomad_conjugation_genes", "rel_coverage", "mobmess_cluster", "mobmess_cluster_type",
            "cluster_members", "cluster_circular_members", "kept_after_mobmess_filter", "replicon_tag",
            "graph_isolated", "graph_connected_to", "graph_repeat_edges", "graph_max_edge_depth_ratio"]
    rows, review = [], []
    for contig, i in info.items():
        m = meta.get(contig, {})
        s = stats.get(m.get("sample", ""), {})
        g = genomad.get(contig, {})
        cl = clusters.get(contig, {})
        mem = members.get(cl.get("cluster"), [])
        row = {
            "contig": contig, "sample": m.get("sample", ""), "dataset": m.get("dataset", ""),
            "length": i["length"], "circular": i["circular"],
            "assembly_contigs": s.get("assembly_contigs", ""),
            "assembly_circular_contigs": s.get("assembly_circular_contigs", ""),
            **rep[contig],
            "rep_genes": operon.get(contig, {}).get("rep_genes", ""),
            "complete_operon": operon.get(contig, {}).get("complete_operon", "0"),
            "genomad_plasmid_score": g.get("plasmid_score", ""), "genomad_hallmarks": g.get("n_hallmarks", ""),
            "genomad_conjugation_genes": g.get("conjugation_genes", ""),
            "rel_coverage": cov.get(contig, {}).get("rel_coverage", ""),
            "mobmess_cluster": cl.get("cluster", ""), "mobmess_cluster_type": cl.get("cluster_type", ""),
            "cluster_members": len(mem), "cluster_circular_members": sum(x["circular"] == "1" for x in mem),
            "kept_after_mobmess_filter": int(contig in kept), "replicon_tag": tags.get(contig, ""),
        }
        gr = graph.get(contig)
        if gr:
            others = [o for o in gr["component_contigs"].split(",") if o != "-"]
            row.update({"graph_isolated": gr["isolated_component"],
                        "graph_connected_to": ",".join(sorted({tags.get(o) or o for o in others})) or "-",
                        "graph_repeat_edges": gr["repeat_edges"],
                        "graph_max_edge_depth_ratio": gr["max_edge_depth_ratio"]})
        else:
            row.update({"graph_isolated": "", "graph_connected_to": "", "graph_repeat_edges": "",
                        "graph_max_edge_depth_ratio": ""})
        rows.append(row)
        if row["circular"] == "1" or row["replicon_tag"] in MAIN:
            continue
        reasons = []
        if row["repA"] + row["repB"] + row["repC"] and contig not in kept:
            reasons.append(
                f"A: linear with RepABC genes ({row['rep_genes'] or 'no repA'}; complete operon: "
                f"{row['complete_operon']}), dropped by MobMess filter ({row['mobmess_cluster_type'] or 'no cluster'}, "
                f"{row['cluster_circular_members']} circular of {row['cluster_members']} members)")
        if gr:
            linked = [tags.get(o, "") for o in gr["component_contigs"].split(",") if o != "-"]
            apart = (not linked) or (tags and not any(t in MAIN for t in linked) and all(linked))
        else:
            apart = False
        n_rep = row["repA"] + row["repB"] + row["repC"]
        markers = ([f"RepABC ({row['rep_genes'] or 'no repA'})"] if n_rep else []) + \
                  ([f"geNomad plasmid (score {row['genomad_plasmid_score']}, {row['genomad_hallmarks']} hallmarks)"]
                   if row["genomad_plasmid_score"] and (int(row["genomad_hallmarks"] or 0) > 0 or
                                                         float(row["genomad_plasmid_score"]) >= args.genomad_min_score)
                   else [])
        if apart and markers:
            where = "alone" if row["graph_isolated"] == "1" else "not linked to a main replicon"
            reasons.append(
                f"B: linear with {' and '.join(markers)}, {where} in the assembly graph "
                f"({row['graph_repeat_edges']} repeat edges, "
                f"max edge depth {row['graph_max_edge_depth_ratio']}x; "
                f"{'kept' if contig in kept else 'dropped'} by MobMess filter)")
        if reasons:
            if row["graph_connected_to"] not in ("", "-"):
                reasons.append(f"graph-connected to {row['graph_connected_to']}")
            review.append(dict(row, review_reason="; ".join(reasons)))

    rows.sort(key=lambda r: (r["sample"], r["contig"]))
    review.sort(key=lambda r: (r["sample"], r["contig"]))
    with open("contig_markers.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open("linear_marker_review.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols + ["review_reason"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(review)
    print(f"{len(rows)} sample contigs; {sum(r['complete_operon'] == '1' for r in rows)} with a complete RepABC operon; "
          f"{len(review)} linear contigs for review")


if __name__ == "__main__":
    main()
