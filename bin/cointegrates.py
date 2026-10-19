#!/usr/bin/env python3
"""
Accessory plasmids that contain another plasmid backbone, and what else they carry.

MobMess lists, for each system, the plasmids that contain its backbone
(<prefix>_systems.txt, compound_plasmids). Many of these backbones are ~3 kb insertion
sequences or NRPS-derived contigs; here only PLASMID partners are kept: the circular
backbone member covering most of the accessory plasmid (MobMess alignment blocks,
<prefix>_ani_blocks.txt) is >= --min-partner-len and the region it covers carries a
replication (Rep) gene. Everything of the accessory plasmid outside its plasmid
partner(s) is an extra segment: candidate integration of another plasmid or MGE.

The publication found the cluster 381 cointegrate (Fig. 5B) by inspecting the MobMess
systems and the Fig. S5A plasmid (a ~73 kb plasmid plus another MGE) with clinker; this
makes both searches systematic.

Per accessory plasmid with a plasmid partner:
  - cointegrate: two plasmid partners from different systems that each cover >= --min-part,
    overlap by <= --max-overlap (of the smaller covered length) and together cover
    >= --min-union of the plasmid. Junctions (backbone overlap and segment ends) are
    reported with the genes within --window bp.
  - extra segments (>= --min-segment bp, not covered by any plasmid partner): genes
    inside (Rep, Par, integrase/recombinase, relaxase/mobilization, conjugation,
    transposase), mobile-element genes within --window bp of each end, and the best
    match among the other MobMess sequences (largest share of the segment aligned;
    sequences of the same MobMess cluster and sequences that contain the whole plasmid are
    skipped): another plasmid carrying the same segment, not necessarily its source. Class, first match: plasmid-derived (Rep inside) > integrative/mobilizable
    (integrase/recombinase or relaxase inside) > transposon-bounded (mobile genes at an end)
    > other.

Outputs: plasmid_partners.tsv, extra_segments.tsv, cointegrates.tsv,
         integration_summary.tsv, cointegrate_junction_genes.tsv, integration_map.png/.svg
"""
import argparse
import csv
import itertools
import os
import re
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)
MOBILE = re.compile(r"transposase|resolvase|integrase|recombinase|insertion element|\bTn3\b", re.I)
REPLICATION = re.compile(r"replication|\brep[ABC]\b", re.I)
PARTITION = re.compile(r"partition|\bpar[AB]\b", re.I)
INTEGRASE = re.compile(r"integrase|recombinase", re.I)
RELAXASE = re.compile(r"relaxase|mobilization", re.I)
CONJUGATION = re.compile(r"type IV secret|conjuga|\bTra[A-Z]\b|\bTrb[A-Z]\b|\bVirB", re.I)
TRANSPOSASE = re.compile(r"transposase|insertion element", re.I)
CLASSES = ["plasmid-derived", "integrative/mobilizable", "transposon-bounded", "other"]
COLORS = {"partner": "#BBBBBB", "plasmid-derived": "#D55E00", "integrative/mobilizable": "#0072B2",
          "transposon-bounded": "#E69F00", "other": "#56B4E9"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clusters", required=True, help="MobMess <prefix>_clusters.txt")
    p.add_argument("--systems", required=True, help="MobMess <prefix>_systems.txt")
    p.add_argument("--blocks", required=True, help="MobMess <prefix>_ani_blocks.txt")
    p.add_argument("--plasmids", required=True, help="accessory plasmid IDs (one per line)")
    p.add_argument("--annotation", required=True, help="Bakta TSV of the accessory plasmids")
    p.add_argument("--min-partner-len", type=int, default=20000, help="plasmid partner: minimum backbone length")
    p.add_argument("--min-part", type=float, default=0.2)
    p.add_argument("--max-overlap", type=float, default=0.2)
    p.add_argument("--min-union", type=float, default=0.8)
    p.add_argument("--min-segment", type=int, default=5000, help="minimum extra segment length (bp)")
    p.add_argument("--window", type=int, default=5000)
    return p.parse_args()


def merge(iv, gap=0):
    out = []
    for a, e in sorted(iv):
        if out and a <= out[-1][1] + gap + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([a, e])
    return out


def length(iv):
    return sum(e - a + 1 for a, e in iv)


def intersect(a, b):
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi >= lo:
            out.append([lo, hi])
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def gaps(iv, xl, circular):
    """Uncovered segments of [1, xl]; on circular plasmids the last gap wraps the origin (end > xl)."""
    iv = merge(iv)
    if not iv:
        return [(1, xl)]
    g = [(iv[k][1] + 1, iv[k + 1][0] - 1) for k in range(len(iv) - 1)]
    if circular:
        g.append((iv[-1][1] + 1, iv[0][0] - 1 + xl))
    else:
        g += [(1, iv[0][0] - 1), (iv[-1][1] + 1, xl)]
    return [x for x in g if x[1] >= x[0]]


def unwrap(a, e, xl):
    """Segment that may cross the origin -> list of plain intervals."""
    return [[a, min(e, xl)]] + ([[1, e - xl]] if e > xl else [])


def is_circular(row):
    return row.get("circular") in ("1", "True", "true")


def main():
    args = parse_args()
    with open(args.plasmids) as fh:
        acc = {l.strip() for l in fh if l.strip()}
    with open(args.clusters) as fh:
        clusters = {r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(args.systems) as fh:
        systems = {r["system_name"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # accessory plasmid -> systems whose backbone it contains (the system file is complete:
    # a plasmid can be the backbone of its own system and a compound of another)
    carried = defaultdict(list)
    for s, r in systems.items():
        for p in (r.get("compound_plasmids") or "").split("|"):
            if p in acc:
                carried[p].append(s)
    X = set(carried)

    # alignment blocks: intervals on each candidate plasmid, per partner sequence
    iv, xlen, plen = defaultdict(list), {}, {}
    with open(args.blocks) as fh:
        r = csv.reader(fh, delimiter="\t")
        h = next(r)
        k = {n: i for i, n in enumerate(h)}
        for row in r:
            ref, qry = row[k["ref"]], row[k["query"]]
            plen[ref], plen[qry] = int(row[k["ref_len"]]), int(row[k["query_len"]])
            if ref == qry:
                continue
            if ref in X:
                a, e = int(row[k["ref_start"]]), int(row[k["ref_end"]])
                iv[(ref, qry)].append((min(a, e), max(a, e)))
                xlen[ref] = int(row[k["ref_len"]])
            if qry in X:
                a, e = int(row[k["query_start"]]), int(row[k["query_end"]])
                iv[(qry, ref)].append((min(a, e), max(a, e)))
                xlen[qry] = int(row[k["query_len"]])
    partners_of = defaultdict(set)
    for x, m in iv:
        partners_of[x].add(m)

    genes = defaultdict(list)
    with open(args.annotation) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and f[0] in X and f[1] == "cds":
                genes[f[0]].append((int(f[2]), int(f[3]), f[5], f[6], f[7]))

    def genes_in(x, intervals):
        """Genes with >= 50% of their length inside the intervals."""
        return [g for g in genes.get(x, [])
                if sum(max(0, min(g[1], e) - max(g[0], a) + 1) for a, e in intervals) >= 0.5 * (g[1] - g[0] + 1)]

    def near(x, points, xl, circular):
        """Genes within --window bp of any point (also across the origin of circular plasmids)."""
        shifts = (0, -xl, xl) if circular else (0,)
        w = args.window
        return [g for g in genes.get(x, [])
                if any(g[1] + sh >= p - w and g[0] + sh <= p + w for p in points for sh in shifts)]

    def count(gs, rx):
        return sum(bool(rx.search(g[4])) for g in gs)

    partner_rows, seg_rows, coint_rows, summary_rows, gene_rows = [], [], [], [], []
    for x in sorted(X):
        xl = xlen.get(x)
        if not xl:
            continue
        circular = is_circular(clusters.get(x, {}))
        # best circular backbone member per system; keep plasmid partners only
        plasmid_partners = []
        for s in carried[x]:
            best = (None, [])
            for m in systems[s]["backbone_plasmids"].split("|"):
                if m == x or not is_circular(clusters.get(m, {})):
                    continue
                cur = merge(iv.get((x, m), []))
                if length(cur) > length(best[1]):
                    best = (m, cur)
            m, cur = best
            if not m or plen.get(m, 0) < args.min_partner_len:
                continue
            if count(genes_in(x, cur), REPLICATION) == 0:
                continue
            plasmid_partners.append((s, m, cur))
        if not plasmid_partners:
            continue
        cl = clusters[x]
        for s, m, cur in plasmid_partners:
            partner_rows.append({
                "plasmid": x, "cluster": cl["cluster"], "cluster_type": cl["cluster_type"], "length": xl,
                "system": s, "backbone_cluster": systems[s]["backbone_cluster"], "partner": m,
                "partner_length": plen[m], "covered_bp": length(cur), "fraction_of_plasmid": round(length(cur) / xl, 3),
                "fraction_of_partner": round(length(cur) / plen[m], 3)})

        # cointegrate: two complementary plasmid partners
        best_pair = None
        for (s1, m1, i1), (s2, m2, i2) in itertools.combinations(plasmid_partners, 2):
            c1, c2 = length(i1) / xl, length(i2) / xl
            ov_iv = intersect(i1, i2)
            ov = length(ov_iv) / max(1, min(length(i1), length(i2)))
            un = length(merge(i1 + i2)) / xl
            if c1 >= args.min_part and c2 >= args.min_part and ov <= args.max_overlap and un >= args.min_union:
                if best_pair is None or un - ov > best_pair[-1]:
                    best_pair = (s1, m1, i1, c1, s2, m2, i2, c2, ov_iv, ov, un, un - ov)
        if best_pair:
            s1, m1, i1, c1, s2, m2, i2, c2, ov_iv, ov, un, _ = best_pair
            ends = [p for a, e in merge(i1, 2000) if e - a >= 5000 for p in (a, e)]
            junctions = merge([list(j) for j in ov_iv] + [[p, p] for p in ends], args.window)
            jg = near(x, [p for a, e in junctions for p in (a, e)], xl, circular)
            coint_rows.append({
                "plasmid": x, "cluster": cl["cluster"], "length": xl,
                "system_1": s1, "backbone_cluster_1": systems[s1]["backbone_cluster"], "backbone_1": m1,
                "covered_kb_1": round(length(i1) / 1000, 1), "fraction_1": round(c1, 3),
                "system_2": s2, "backbone_cluster_2": systems[s2]["backbone_cluster"], "backbone_2": m2,
                "covered_kb_2": round(length(i2) / 1000, 1), "fraction_2": round(c2, 3),
                "overlap": round(ov, 3), "union": round(un, 3),
                "junctions": ";".join(f"{a}-{e}" for a, e in junctions),
                "mobile_genes_at_junctions": ";".join(f"{g[2]}:{g[4]}" for g in jg if MOBILE.search(g[4]))})
            gene_rows += [{"type": "cointegrate junction", "plasmid": x, "start": g[0], "end": g[1], "locus": g[2],
                           "gene": g[3], "product": g[4], "mobile": int(bool(MOBILE.search(g[4])))} for g in jg]

        # extra segments: not covered by any plasmid partner
        core = merge([list(t) for _, _, cur in plasmid_partners for t in cur])
        own_cluster = {c for c, r in clusters.items() if r["cluster"] == cl["cluster"]}
        partner_seqs = {m for _, m, _ in plasmid_partners}
        segs = [g for g in gaps(core, xl, circular) if g[1] - g[0] + 1 >= args.min_segment]
        classes = []
        for a, e in segs:
            seg = unwrap(a, e, xl)
            e_pos = e if e <= xl else e - xl
            inside = genes_in(x, seg)
            flank = near(x, [a, e_pos], xl, circular)
            mobile_flank = [g for g in flank if MOBILE.search(g[4]) and g not in inside]
            # best match: other sequence aligning to most of the segment (not containing the whole plasmid)
            donor, dcov = None, 0
            for m in partners_of[x]:
                if m in own_cluster or m in partner_seqs:
                    continue
                mi = merge(iv[(x, m)])
                if length(mi) >= 0.8 * xl:
                    continue
                c = length(intersect(mi, merge(seg)))
                if c > dcov:
                    donor, dcov = m, c
            n_rep, n_int, n_rel = count(inside, REPLICATION), count(inside, INTEGRASE), count(inside, RELAXASE)
            cls = ("plasmid-derived" if n_rep else "integrative/mobilizable" if (n_int or n_rel)
                   else "transposon-bounded" if mobile_flank else "other")
            classes.append(cls)
            seg_rows.append({
                "plasmid": x, "cluster": cl["cluster"], "length": xl, "start": a, "end": e_pos, "segment_bp": e - a + 1,
                "class": cls, "genes": len(inside), "replication": n_rep, "partition": count(inside, PARTITION),
                "integrase_recombinase": n_int, "relaxase_mobilization": n_rel,
                "conjugation": count(inside, CONJUGATION), "transposase": count(inside, TRANSPOSASE),
                "mobile_genes_at_ends": ";".join(f"{g[2]}:{g[4]}" for g in mobile_flank),
                "best_match": donor or "", "best_match_cluster": clusters.get(donor, {}).get("cluster", ""),
                "best_match_type": clusters.get(donor, {}).get("cluster_type", ""),
                "best_match_fraction_of_segment": round(dcov / (e - a + 1), 3)})
            gene_rows += [{"type": f"segment end ({cls})", "plasmid": x, "start": g[0], "end": g[1], "locus": g[2],
                           "gene": g[3], "product": g[4], "mobile": 1} for g in mobile_flank]

        summary_rows.append({
            "plasmid": x, "cluster": cl["cluster"], "cluster_type": cl["cluster_type"], "length": xl,
            "plasmid_partners": len(plasmid_partners),
            "partner_systems": "|".join(s for s, _, _ in plasmid_partners),
            "partner_backbone_clusters": "|".join(systems[s]["backbone_cluster"] for s, _, _ in plasmid_partners),
            "covered_by_partners_bp": length(core), "extra_segments": len(segs),
            "extra_bp": sum(e - a + 1 for a, e in segs), "cointegrate": int(bool(best_pair)),
            **{f"n_{c.split('/')[0].replace('-', '_')}": classes.count(c) for c in CLASSES}})

    def write(path, rows, cols):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    write("plasmid_partners.tsv", partner_rows,
          ["plasmid", "cluster", "cluster_type", "length", "system", "backbone_cluster", "partner", "partner_length",
           "covered_bp", "fraction_of_plasmid", "fraction_of_partner"])
    write("extra_segments.tsv", seg_rows,
          ["plasmid", "cluster", "length", "start", "end", "segment_bp", "class", "genes", "replication", "partition",
           "integrase_recombinase", "relaxase_mobilization", "conjugation", "transposase", "mobile_genes_at_ends",
           "best_match", "best_match_cluster", "best_match_type", "best_match_fraction_of_segment"])
    write("cointegrates.tsv", coint_rows,
          ["plasmid", "cluster", "length", "system_1", "backbone_cluster_1", "backbone_1", "covered_kb_1", "fraction_1",
           "system_2", "backbone_cluster_2", "backbone_2", "covered_kb_2", "fraction_2", "overlap", "union",
           "junctions", "mobile_genes_at_junctions"])
    write("integration_summary.tsv", summary_rows,
          ["plasmid", "cluster", "cluster_type", "length", "plasmid_partners", "partner_systems",
           "partner_backbone_clusters", "covered_by_partners_bp", "extra_segments", "extra_bp", "cointegrate"]
          + [f"n_{c.split('/')[0].replace('-', '_')}" for c in CLASSES])
    write("cointegrate_junction_genes.tsv", gene_rows, ["type", "plasmid", "start", "end", "locus", "gene", "product", "mobile"])

    plot(summary_rows, partner_rows, seg_rows)
    n_cl = len({r["cluster"] for r in summary_rows})
    print(f"{len(summary_rows)} accessory plasmids ({n_cl} clusters) contain a plasmid partner backbone; "
          f"{len(seg_rows)} extra segments >= {args.min_segment} bp; cointegrates: {len(coint_rows)} plasmids in "
          f"{len({r['cluster'] for r in coint_rows})} clusters")


def plot(summary_rows, partner_rows, seg_rows):
    """One bar per plasmid (sorted by length): plasmid-partner coverage and extra segments by class."""
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    rows = sorted(summary_rows, key=lambda r: r["length"])
    fig, ax = plt.subplots(figsize=(7, max(3, 0.06 * len(rows) + 1)))
    if rows:
        segs = defaultdict(list)
        for s in seg_rows:
            segs[s["plasmid"]].append(s)
        for i, r in enumerate(rows):
            ax.barh(i, r["length"] / 1000, color=COLORS["partner"], height=0.8, linewidth=0)
            for s in segs[r["plasmid"]]:
                for a, e in unwrap(s["start"], s["start"] + s["segment_bp"] - 1, r["length"]):
                    ax.barh(i, (e - a + 1) / 1000, left=(a - 1) / 1000, color=COLORS[s["class"]], height=0.8, linewidth=0)
            if r["cointegrate"]:
                ax.text(r["length"] / 1000 + 3, i, "cointegrate", va="center", fontsize=5)
    ax.set_yticks([])
    ax.set_ylabel(f"Accessory plasmids with a plasmid partner (n = {len(rows)})")
    ax.set_xlabel("Position (kb)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=[Patch(color=COLORS["partner"], label="covered by plasmid partner(s)")] +
              [Patch(color=COLORS[c], label=f"extra: {c}") for c in CLASSES], frameon=False, fontsize=6, loc="lower right")
    fig.savefig("integration_map.png", dpi=300, bbox_inches="tight")
    fig.savefig("integration_map.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
