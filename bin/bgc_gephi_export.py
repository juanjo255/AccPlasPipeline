#!/usr/bin/env python3
"""
Export the BiG-SCAPE mix-bin network to Gephi node/edge CSVs, annotated with BGC novelty
(bin/bgc_novelty.py) and the plasmid each BGC comes from (final MobMess cluster, replicon group).

Outputs:
  nodes.csv    Id, Label, source (query/MIBiG), antiSMASH class, BiG-SCAPE family, connected component,
               plasmid, plasmid_cluster (all + top-N, others 'other'), replicon_group, known/novel label,
               BiG-SLiCE GCF with its novelty / size / match, per-BGC antiSMASH-DB distance and match
  edges.csv    Source, Target, Type, Weight (= 1 - distance), distance, jaccard, adjacency, dss, bridge
               (bridge = edge joining two different BiG-SCAPE families)
  bridges.tsv  family pairs joined by bridge edges, strongest first
Import in Gephi: File > Import spreadsheet (nodes first, then edges, append to the same workspace).
"""
import argparse
import csv
import re
from collections import defaultdict

MIBIG_RE = re.compile(r"^BGC\d{7}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--network", required=True, help="BiG-SCAPE mix_c<cutoff>.network")
    p.add_argument("--clustering", required=True, help="BiG-SCAPE mix_clustering_c<cutoff>.tsv")
    p.add_argument("--annotations", required=True, help="BiG-SCAPE record_annotations.tsv")
    p.add_argument("--bgc-novelty", required=True, help="bgc_novelty.tsv (bin/bgc_novelty.py)")
    p.add_argument("--gcf-novelty", required=True, help="gcf_novelty_bigslice.tsv (bin/bgc_novelty.py)")
    p.add_argument("--plasmid-clusters", required=True, help="MobMess <prefix>_clusters.txt (contig, cluster)")
    p.add_argument("--replicon-groups", default=None, help="repA.replicon_groups.tsv (sequence, replicon_group)")
    p.add_argument("--top-clusters", type=int, default=8, help="plasmid clusters with most BGCs kept by name")
    return p.parse_args()


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    a = parse_args()
    ann = {r["GBK"]: r for r in read_tsv(a.annotations)}
    fam = {r["GBK"]: r for r in read_tsv(a.clustering)}
    edges = read_tsv(a.network)
    labels = {r["bgc"]: r for r in read_tsv(a.bgc_novelty)}
    gcfs = {r["bigslice_gcf"]: r for r in read_tsv(a.gcf_novelty)}
    size_col = next((c for c in (read_tsv(a.gcf_novelty)[:1] or [{}])[0] if c.startswith("size_ge")), "size_ge5")
    pls_cluster = {r["contig"]: r["cluster"] for r in read_tsv(a.plasmid_clusters)}
    rep_groups = defaultdict(set)
    if a.replicon_groups:
        for r in read_tsv(a.replicon_groups):
            if r.get("replicon_group"):
                rep_groups[r["sequence"]].add(r["replicon_group"])

    # Nodes = every record in the mix-bin network (query BGCs + MIBiG BGCs linked to them)
    ids = set(fam) | {e["GBK_a"] for e in edges} | {e["GBK_b"] for e in edges}
    # Plasmid of origin = BGC id without ".regionNNN"; largest plasmid clusters keep their id
    record_of = {k: re.sub(r"\.region\d+$", "", k) for k in ids if not MIBIG_RE.match(k)}
    cl_counts = defaultdict(int)
    for rec in record_of.values():
        if rec in pls_cluster:
            cl_counts[pls_cluster[rec]] += 1
    top = {c for c, _ in sorted(cl_counts.items(), key=lambda x: -x[1])[:a.top_clusters]}

    node_rows = []
    for k in sorted(ids):
        an, f, lab = ann.get(k, {}), fam.get(k, {}), labels.get(k, {})
        g = gcfs.get(lab.get("bigslice_gcf", ""), {})
        is_mibig = bool(MIBIG_RE.match(k))
        rec = record_of.get(k, "")
        cl = pls_cluster.get(rec, "")
        node_rows.append({
            "Id": k,
            # MIBiG nodes get their description (usually names the compound) as label
            "Label": f"{k} | {an.get('Description', '')[:60]}" if is_mibig else k,
            "source": "MIBiG" if is_mibig else "query",
            "antismash_class": an.get("Class", ""), "category": an.get("Category", ""),
            "organism": an.get("Organism", ""),
            "plasmid": rec,
            # prefixed so Gephi imports the MobMess cluster as a category, not a number
            "plasmid_cluster": f"C{cl}" if cl else "",
            "plasmid_cluster_top": (f"C{cl}" if cl in top else "other") if rec else "",
            "replicon_group": ";".join(sorted(rep_groups.get(rec, ()))),
            "bigscape_family": f.get("Family", ""), "connected_component": f.get("CC", ""),
            "novelty_label": "MIBiG_reference" if is_mibig else lab.get("label", ""),
            "bigslice_gcf": lab.get("bigslice_gcf", ""),
            "gcf_novelty": g.get("novelty_mean_asdb4_dist", ""), "gcf_size": g.get("n_bgcs", ""),
            "gcf_size_ge5": g.get(size_col, ""),
            "gcf_antismashdb_match": g.get("antismashdb_match", ""),
            "gcf_frac_antismashdb_match": g.get("frac_members_antismashdb_match", ""),
            "antismashdb_match": lab.get("antismashdb_match", ""),
            "asdb4_best_dist": lab.get("asdb4_best_dist", ""),
            "on_contig_edge": lab.get("on_contig_edge", ""),
            "nearest_mibig": lab.get("nearest_mibig", ""), "nearest_mibig_dist": lab.get("nearest_mibig_dist", ""),
        })

    edge_rows, bridges = [], defaultdict(list)
    for e in edges:
        d = float(e["distance"])
        fa, fb = fam.get(e["GBK_a"], {}).get("Family", ""), fam.get(e["GBK_b"], {}).get("Family", "")
        bridge = bool(fa and fb and fa != fb)
        edge_rows.append({"Source": e["GBK_a"], "Target": e["GBK_b"], "Type": "Undirected",
                          "Weight": round(1 - d, 4), "distance": d, "jaccard": e["jaccard"],
                          "adjacency": e["adjacency"], "dss": e["dss"], "bridge": bridge})
        if bridge:
            bridges[tuple(sorted((fa, fb)))].append((d, e["GBK_a"], e["GBK_b"], e["jaccard"], e["adjacency"], e["dss"]))

    for name, rows in (("nodes.csv", node_rows), ("edges.csv", edge_rows)):
        with open(name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    # Family pairs joined by at least one edge below the cutoff (the long edges in the viewer)
    cls = {r["Id"]: r["antismash_class"] for r in node_rows}
    with open("bridges.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family_a", "family_b", "n_edges", "min_distance", "best_pair_a", "class_a",
                    "best_pair_b", "class_b", "jaccard", "adjacency", "dss"])
        for (fa, fb), lst in sorted(bridges.items(), key=lambda x: min(x[1])[0]):
            d, ga, gb, jac, adj, dss = min(lst)
            w.writerow([fa, fb, len(lst), d, ga, cls.get(ga, ""), gb, cls.get(gb, ""), jac, adj, dss])

    n_q = sum(r["source"] == "query" for r in node_rows)
    print(f"{len(node_rows)} nodes ({n_q} query, {len(node_rows) - n_q} MIBiG), {len(edge_rows)} edges, "
          f"{sum(r['bridge'] for r in edge_rows)} bridge edges; "
          f"{sum(bool(r['plasmid_cluster']) for r in node_rows if r['source'] == 'query')}/{n_q} BGCs with a plasmid cluster")


if __name__ == "__main__":
    main()
