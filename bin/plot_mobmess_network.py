#!/usr/bin/env python3
"""
Interactive MobMess containment network (all clusters, before filtering).

Nodes are sequences (MobMess <prefix>_contig_graph.graphml), edges are
containment relationships. Nodes are coloured by MobMess cluster type and sized
by sequence length, as in the publication (Fig. S2A).

Node attributes added to the MobMess graph:
  length, gc         from the MobMess input FASTA (GC = (G+C) / length * 100)
  source             sample / plsdb
  genomad            plasmid / not_plasmid (sample contigs), not_tested (PLSDB)
  passed_filter      1 if the sequence passed the MobMess cluster filter, else 0

Outputs:
  <prefix>.html      interactive network (ipysigma; open in a web browser)
  <prefix>.graphml   annotated graph (e.g. for Gephi or Cytoscape)
  <prefix>.nodes.tsv node attribute table
"""
import argparse
import csv
import sys

import networkx as nx
from ipysigma import Sigma

# Same colours as the publication notebook
CLUSTER_TYPE_COLORS = {
    "backbone":              "#2196F3",
    "maximal_not_in_system": "#FF9800",
    "fragment":              "#E91E8C",
    "compound":              "#00BCD4",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--graph", required=True, help="MobMess <prefix>_contig_graph.graphml")
    p.add_argument("--fasta", required=True, help="MobMess input FASTA")
    p.add_argument("--plsdb-ids", required=True, help="PLSDB accessions in the MobMess input")
    p.add_argument("--filtered-ids", required=True, help="sequences that passed the cluster filter")
    p.add_argument("--genomad", nargs="*", default=[], help="geNomad <sample>_plasmid_summary.tsv files")
    p.add_argument("--prefix", default="mobmess_network")
    return p.parse_args()


def length_gc(path):
    # Intent: length and GC% per sequence ID (first word of the header)
    stats, seq_id, n, gc = {}, None, 0, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    stats[seq_id] = (n, gc)
                seq_id, n, gc = line[1:].split(None, 1)[0], 0, 0
            else:
                s = line.strip().upper()
                n += len(s)
                gc += s.count("G") + s.count("C")
    if seq_id is not None:
        stats[seq_id] = (n, gc)
    return {k: (n, round(gc / n * 100, 2) if n else 0.0) for k, (n, gc) in stats.items()}


def read_ids(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def main():
    args = parse_args()
    graph = nx.read_graphml(args.graph)
    stats = length_gc(args.fasta)
    plsdb = read_ids(args.plsdb_ids)
    passed = read_ids(args.filtered_ids)

    genomad_plasmids = set()
    for path in args.genomad:
        with open(path, newline="") as fh:
            genomad_plasmids.update(row["seq_name"] for row in csv.DictReader(fh, delimiter="\t"))

    names = nx.get_node_attributes(graph, "name")
    missing = [n for n in names.values() if n not in stats]
    if missing:
        sys.exit(f"ERROR: graph nodes not in {args.fasta}: {', '.join(missing[:5])}")
    if len(names) != len(stats):
        sys.exit(f"ERROR: graph has {len(names)} nodes but {args.fasta} has {len(stats)} sequences")

    for node, name in names.items():
        data = graph.nodes[node]
        is_plsdb = name in plsdb
        data["cluster"] = int(data["cluster"])
        data["circular"] = int(data["circular"])
        data["length"], data["gc"] = stats[name]
        data["source"] = "plsdb" if is_plsdb else "sample"
        data["genomad"] = "not_tested" if is_plsdb else ("plasmid" if name in genomad_plasmids else "not_plasmid")
        data["passed_filter"] = int(name in passed)

    # Interactive HTML; the ForceAtlas2 layout runs in the browser
    Sigma.write_html(
        graph,
        f"{args.prefix}.html",
        fullscreen=True,
        node_color="cluster_type",
        node_color_palette=CLUSTER_TYPE_COLORS,
        node_size="length",
        node_size_range=(3, 20),
        node_label="name",
        node_border_color_from="node",
        default_edge_type="curve",
        layout_settings={"scalingRatio": 0.8},
        label_density=0,
    )

    nx.write_graphml(graph, f"{args.prefix}.graphml")

    cols = ["name", "cluster", "cluster_type", "circular", "systems", "length", "gc",
            "source", "genomad", "passed_filter"]
    with open(f"{args.prefix}.nodes.tsv", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(cols)
        for node in sorted(graph.nodes, key=lambda n: graph.nodes[n]["name"]):
            writer.writerow([graph.nodes[node].get(c, "") for c in cols])

    print(f"Network: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges, "
          f"{sum(graph.nodes[n]['passed_filter'] for n in graph)} passed the filter")


if __name__ == "__main__":
    main()
