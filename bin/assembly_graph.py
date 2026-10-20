#!/usr/bin/env python3
"""
Assembly-graph context of each contig of a Flye assembly.

Connected components of the Flye graph (assembly_graph.gfa, S/L lines) are mapped to the
contigs through the graph_path column of assembly_info.txt. A linear contig that is alone in
its component, at a depth close to the chromosome, is probably a complete replicon that Flye
could not circularize, typically because of repeats (IS elements, transposons) that the path
traverses several times. A contig connected to other contigs shares graph sequence with them
(for the chromosome: likely a fragment or an integrated element; for pSymA/pSymB: ambiguous).

Repeat edges of a contig: edges its path uses more than once, or with depth >= --repeat-ratio
times the contig depth (Flye 'cov.').

Output <sample>.graph.tsv: contig, graph_component, component_contigs (other contigs in the
component), isolated_component, component_bp, component_edges, path_edges, repeat_edges,
max_edge_depth_ratio
"""
import argparse
from collections import Counter, defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gfa", required=True, help="Flye assembly_graph.gfa")
    p.add_argument("--assembly-info", required=True, help="Flye assembly_info.txt")
    p.add_argument("--sample", required=True)
    p.add_argument("--add-sample-suffix", action="store_true",
                   help="append _<sample> to contig names unless already present (as Step 1)")
    p.add_argument("--repeat-ratio", type=float, default=1.8)
    return p.parse_args()


def main():
    args = parse_args()
    length, depth, adj = {}, {}, defaultdict(set)
    with open(args.gfa) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                tags = dict(t.split(":", 2)[::2] for t in f[3:] if t.count(":") >= 2)
                length[f[1]] = len(f[2]) if f[2] != "*" else int(tags.get("LN", 0))
                depth[f[1]] = float(tags.get("dp", 0))
            elif f[0] == "L":
                adj[f[1]].add(f[3])
                adj[f[3]].add(f[1])

    comp, n = {}, 0
    for s in length:
        if s in comp:
            continue
        n += 1
        stack = [s]
        while stack:
            x = stack.pop()
            if x not in comp:
                comp[x] = n
                stack.extend(adj[x] - comp.keys())
    comp_edges = defaultdict(list)
    for s, c in comp.items():
        comp_edges[c].append(s)

    contigs = []
    with open(args.assembly_info) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            name = f[0]
            if args.add_sample_suffix and not name.endswith("_" + args.sample):
                name = f"{name}_{args.sample}"
            # graph_path: comma-separated edge numbers with optional sign; '*' marks gaps
            edges = [f"edge_{e.strip('+-')}" for e in f[7].split(",") if e.strip("+-*")]
            contigs.append((name, float(f[2]), edges))

    by_comp = defaultdict(list)
    for name, _, edges in contigs:
        for c in {comp[e] for e in edges if e in comp}:
            by_comp[c].append(name)

    with open(f"{args.sample}.graph.tsv", "w") as out:
        out.write("contig\tgraph_component\tcomponent_contigs\tisolated_component\tcomponent_bp\t"
                  "component_edges\tpath_edges\trepeat_edges\tmax_edge_depth_ratio\n")
        for name, cov, edges in contigs:
            cs = sorted({comp[e] for e in edges if e in comp})
            others = sorted({o for c in cs for o in by_comp[c] if o != name})
            used = Counter(edges)
            rep = {e for e in used if used[e] > 1 or (cov and depth.get(e, 0) >= args.repeat_ratio * cov)}
            ratio = max((depth.get(e, 0) / cov for e in used), default=0) if cov else 0
            out.write("\t".join(map(str, [
                name, ",".join(map(str, cs)), ",".join(others) or "-", int(not others),
                sum(length[e] for c in cs for e in comp_edges[c]), sum(len(comp_edges[c]) for c in cs),
                len(used), len(rep), f"{ratio:.1f}"])) + "\n")


if __name__ == "__main__":
    main()
