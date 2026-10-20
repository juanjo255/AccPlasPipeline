#!/usr/bin/env python3
"""
Compound plasmids: what do they carry that their parent plasmid does not?

For every accessory plasmid in a MobMess *compound* cluster (Step 4 run, whose cluster IDs are
the ones in the publication), the parents are the plasmid backbones of its systems: for each
system, the circular backbone member covering most of the compound (MobMess alignment blocks),
kept as a parent when it is >= --min-partner-len and its covered region carries a replication
gene. Backbones that are not plasmids (~3 kb IS circles, NRPS-derived contigs) are reported but
not used as parents.

Cointegrate: two parents from different systems that each cover >= --min-part of the compound,
overlap by <= --max-overlap (of the smaller covered length) and together cover >= --min-union.

Extra sequence: every region >= --min-segment bp not covered by any parent (circular-aware). Its
content is read from several sources, all reported as hints (ICE / IME / phage annotations are
not well established):
  replication    RepABC HMM genes (Step 4b), MOB-typer replicon hits, Bakta replication genes
  relaxase       MOB-typer relaxase, geNomad CONJscan MOB*, Bakta relaxase / mobilization names
  MPF / T4SS     MOB-typer mate-pair-formation, geNomad CONJscan T4SS genes, Bakta conjugation names
  integrase      Bakta integrase / site-specific recombinase names (not Tn3 resolvases)
  transposase    Bakta transposase / IS names
  phage          geNomad proviruses overlapping the region, geNomad virus hallmark genes, Bakta phage
                 structural gene names (terminase, capsid, portal, tail, baseplate, holin ...)
Class (first match):
  phage-plasmid               replication + phage
  plasmid-derived (second parent absent)   replication, no other sequence covers >= 80% of the region:
                              the region likely came from a second plasmid that is not in the dataset
                              (the first parent, the backbone, always is)
  plasmid-derived (second parent present)  replication, another sequence covers >= 80% of the region
  phage                       provirus covering >= 50% of the provirus, or >= 2 phage structural genes
  ICE-like                    integrase + relaxase + MPF
  IME-like                    integrase + relaxase
  integrative element         integrase inside or at an end, no relaxase found
  MOB/MPF module              relaxase and/or MPF without integrase
  transposon acquisition      transposase inside or at an end
  other
The best match among the other MobMess sequences (largest share of the region aligned; same
cluster and sequences containing the whole plasmid skipped) shows whether the extra sequence
exists elsewhere in the data.

Outputs: compound_parents.tsv, extra_segments.tsv, cointegrates.tsv, integration_summary.tsv,
         cointegrate_junction_genes.tsv, integration_map.png/.svg
"""
import argparse
import csv
import itertools
import os
import re
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)
REPLICATION = re.compile(r"replication|\brep[ABC]\b", re.I)
RELAXASE = re.compile(r"relaxase|mobilization|\bmob[A-Z]\b", re.I)
MPF = re.compile(r"type IV secret|conjugal|conjugative transfer|\bTra[A-Z]\b|\bTrb[A-Z]\b|\bVirB\d*", re.I)
INTEGRASE = re.compile(r"integrase|site-specific (DNA )?recombinase|recombinase", re.I)
RESOLVASE = re.compile(r"resolvase|Tn3", re.I)
TRANSPOSASE = re.compile(r"transposase|insertion (sequence|element)|\bIS\d", re.I)
# Phage structural / virion genes only ('phage shock protein', 'anti-phage' or 'phage integrase' are not phages)
PHAGE = re.compile(r"terminase|capsid|portal protein|tail (fiber|fibre|tube|sheath|tape|assembly|protein)|"
                   r"baseplate|prohead|head-tail|head (morphogenesis|decoration|completion)|holin|virion", re.I)
MOBILE = re.compile(r"transposase|resolvase|integrase|recombinase|insertion element|\bTn3\b", re.I)
CLASSES = ["phage-plasmid", "plasmid-derived (second parent absent)", "plasmid-derived (second parent present)", "phage", "ICE-like", "IME-like",
           "integrative element", "MOB/MPF module", "transposon acquisition", "other"]
COLORS = {"parent": "#BBBBBB", "phage-plasmid": "#8c564b", "plasmid-derived (second parent absent)": "#D55E00",
          "plasmid-derived (second parent present)": "#f0a070", "phage": "#9467bd", "ICE-like": "#0072B2", "IME-like": "#56B4E9", "integrative element": "#009E73",
          "MOB/MPF module": "#CC79A7", "transposon acquisition": "#E69F00", "other": "#e0e0e0"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clusters", required=True, help="MobMess <prefix>_clusters.txt")
    p.add_argument("--systems", required=True, help="MobMess <prefix>_systems.txt")
    p.add_argument("--blocks", required=True, help="MobMess <prefix>_ani_blocks.txt")
    p.add_argument("--plasmids", required=True, help="accessory plasmid IDs (one per line)")
    p.add_argument("--annotation", required=True, help="Bakta TSV of the accessory plasmids")
    p.add_argument("--repabc", nargs="*", default=[], help="RepABC HMM genes (Step 4b, Bakta layout)")
    p.add_argument("--mob-biomarker", nargs="*", default=[], help="MOB-typer *_mobtyper_biomarker.tsv")
    p.add_argument("--genomad-genes", nargs="*", default=[], help="geNomad *_plasmid_genes.tsv")
    p.add_argument("--genomad-virus", nargs="*", default=[], help="geNomad *_virus_summary.tsv (proviruses)")
    p.add_argument("--min-partner-len", type=int, default=20000, help="parent: minimum backbone length")
    p.add_argument("--min-part", type=float, default=0.2)
    p.add_argument("--max-overlap", type=float, default=0.2)
    p.add_argument("--min-union", type=float, default=0.8)
    p.add_argument("--min-segment", type=int, default=5000, help="minimum extra region (bp)")
    p.add_argument("--window", type=int, default=5000, help="genes within this distance of a region end")
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
    return [[a, min(e, xl)]] + ([[1, e - xl]] if e > xl else [])


def inside(s, t, region, frac=0.5):
    """Feature [s, t] has >= frac of its length inside the region intervals."""
    return sum(max(0, min(t, e) - max(s, a) + 1) for a, e in region) >= frac * (t - s + 1)


def is_circular(row):
    return row.get("circular") in ("1", "True", "true")


def read_tsv(paths):
    for path in paths:
        with open(path, newline="") as fh:
            yield from csv.DictReader((l for l in fh if not l.startswith("#") or l.startswith("#Sequence")),
                                      delimiter="\t")


def main():
    args = parse_args()
    with open(args.plasmids) as fh:
        acc = {l.strip() for l in fh if l.strip()}
    with open(args.clusters) as fh:
        clusters = {r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(args.systems) as fh:
        systems = {r["system_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    compounds = [c for c, r in clusters.items() if c in acc and r["cluster_type"] == "compound" and r.get("systems")]
    X = set(compounds)

    # MobMess alignment blocks: intervals on each compound, per partner sequence
    iv, xlen, plen = defaultdict(list), {}, {}
    with open(args.blocks) as fh:
        r = csv.reader(fh, delimiter="\t")
        k = {n: i for i, n in enumerate(next(r))}
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

    # Evidence per compound: (start, end, label) lists
    genes = defaultdict(list)   # Bakta: (start, end, locus, gene, product)
    with open(args.annotation) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and f[0] in X and f[1] == "cds":
                genes[f[0]].append((int(f[2]), int(f[3]), f[5], f[6], f[7]))
    ev = defaultdict(lambda: defaultdict(list))   # contig -> kind -> [(start, end, label)]
    for r in read_tsv(args.repabc):
        if r["#Sequence Id"] in X:
            ev[r["#Sequence Id"]]["replication"].append((int(r["Start"]), int(r["Stop"]), f"{r['Gene']} (HMM)"))
    for r in read_tsv(args.mob_biomarker):
        c = r.get("sseqid")
        if c in X:
            s, e = sorted((int(r["sstart"]), int(r["send"])))
            kind = {"replicon": "replication", "relaxase": "relaxase", "mate-pair-formation": "MPF",
                    "oriT": "relaxase"}.get(r["biomarker"])
            if kind:
                ev[c][kind].append((s, e, f"{r['biomarker']} (MOB-typer)"))
    for r in read_tsv(args.genomad_genes):
        c = r["gene"].rsplit("_", 1)[0]
        cs = r.get("annotation_conjscan", "NA")
        if c in X and cs not in ("NA", ""):
            kind = "relaxase" if cs.upper().startswith("MOB") else "MPF"
            ev[c][kind].append((int(r["start"]), int(r["end"]), f"{cs} (geNomad)"))
        if c in X and r.get("virus_hallmark") == "1":
            ev[c]["phage gene"].append((int(r["start"]), int(r["end"]),
                                        f"{(r.get('annotation_description') or 'virus hallmark')[:40]} (geNomad hallmark)"))
    for r in read_tsv(args.genomad_virus):
        c, _, tag = r["seq_name"].partition("|provirus_")
        if c in X and tag:
            s, e = map(int, tag.split("_")[:2])
            ev[c]["provirus"].append((s, e, f"provirus {s}-{e} (geNomad, score {r.get('virus_score', '?')})"))
    for c, gs in genes.items():
        for s, e, locus, g, prod in gs:
            if REPLICATION.search(prod):
                ev[c]["replication"].append((s, e, f"{g or prod} (Bakta)"))
            if RELAXASE.search(prod):
                ev[c]["relaxase"].append((s, e, f"{g or prod} (Bakta)"))
            if MPF.search(prod):
                ev[c]["MPF"].append((s, e, f"{g or prod} (Bakta)"))
            if INTEGRASE.search(prod) and not RESOLVASE.search(prod):
                ev[c]["integrase"].append((s, e, f"{g or prod} (Bakta)"))
            if TRANSPOSASE.search(prod) or RESOLVASE.search(prod):
                ev[c]["transposase"].append((s, e, f"{g or prod} (Bakta)"))
            if PHAGE.search(prod):
                ev[c]["phage gene"].append((s, e, f"{g or prod} (Bakta)"))

    def best_member(x, system):
        best = (None, [])
        for m in systems[system]["backbone_plasmids"].split("|"):
            if m == x or not is_circular(clusters.get(m, {})):
                continue
            cur = merge(iv.get((x, m), []))
            if length(cur) > length(best[1]):
                best = (m, cur)
        return best

    def near(x, points, xl, circular):
        shifts = (0, -xl, xl) if circular else (0,)
        w = args.window
        return [g for g in genes.get(x, [])
                if any(g[1] + sh >= p - w and g[0] + sh <= p + w for p in points for sh in shifts)]

    parent_rows, seg_rows, coint_rows, summary_rows, gene_rows = [], [], [], [], []
    for x in sorted(compounds):
        xl = xlen.get(x)
        if not xl:
            continue
        cl = clusters[x]
        circular = is_circular(cl)
        parents, others = [], []
        for s in cl["systems"].split("|"):
            if s not in systems:
                continue
            m, cur = best_member(x, s)
            if not m:
                continue
            has_rep = any(inside(a, e, cur) for a, e, _ in ev[x]["replication"])
            if plen.get(m, 0) >= args.min_partner_len and has_rep:
                parents.append((s, m, cur))
            else:
                others.append(f"{s}:{systems[s]['backbone_cluster']} ({plen.get(m, 0) // 1000} kb)")
        parent_rows.append({
            "plasmid": x, "cluster": cl["cluster"], "length": xl, "circular": int(circular), "systems": cl["systems"],
            "plasmid_parents": "|".join(f"{s}:{systems[s]['backbone_cluster']}:{m} ({length(c)} bp)" for s, m, c in parents),
            "non_plasmid_backbones": "|".join(others),
            "note": "" if parents else "no plasmid parent (only small / non-plasmid backbones)"})
        if not parents:
            continue

        # Cointegrate: two complementary parents
        best_pair = None
        for (s1, m1, i1), (s2, m2, i2) in itertools.combinations(parents, 2):
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
                           "gene": g[3], "product": g[4]} for g in jg if MOBILE.search(g[4])]

        # Extra sequence: not covered by any parent
        core = merge([list(t) for _, _, cur in parents for t in cur])
        own = {c for c, r in clusters.items() if r["cluster"] == cl["cluster"]}
        parent_seqs = {m for _, m, _ in parents}
        classes = []
        for a, e in (g for g in gaps(core, xl, circular) if g[1] - g[0] + 1 >= args.min_segment):
            region = unwrap(a, e, xl)
            e_pos = e if e <= xl else e - xl
            seg_bp = e - a + 1
            ends = near(x, [a, e_pos], xl, circular)
            hit = {kind: [lab for s_, t_, lab in ev[x][kind] if inside(s_, t_, region)]
                   for kind in ("replication", "relaxase", "MPF", "integrase", "transposase", "phage gene")}
            pro = [(s_, t_, lab) for s_, t_, lab in ev[x]["provirus"]
                   if sum(max(0, min(t_, re_) - max(s_, ra) + 1) for ra, re_ in region) >= 0.5 * (t_ - s_ + 1)]
            int_end = [f"{g[3] or g[4]} (Bakta, end)" for g in ends
                       if INTEGRASE.search(g[4]) and not RESOLVASE.search(g[4]) and not inside(g[0], g[1], region)]
            tn_end = [f"{g[3] or g[4]} (Bakta, end)" for g in ends
                      if (TRANSPOSASE.search(g[4]) or RESOLVASE.search(g[4])) and not inside(g[0], g[1], region)]
            donor, dcov = None, 0
            for m in partners_of[x]:
                if m in own or m in parent_seqs:
                    continue
                mi = merge(iv[(x, m)])
                if length(mi) >= 0.8 * xl:
                    continue
                c = length(intersect(mi, merge(region)))
                if c > dcov:
                    donor, dcov = m, c
            in_data = dcov >= 0.8 * seg_bp
            phage = bool(pro) or len(hit["phage gene"]) >= 2
            integ = bool(hit["integrase"] or int_end)
            if hit["replication"] and phage:
                cls = "phage-plasmid"
            elif hit["replication"]:
                cls = "plasmid-derived (second parent present)" if in_data else "plasmid-derived (second parent absent)"
            elif phage:
                cls = "phage"
            elif integ and hit["relaxase"] and hit["MPF"]:
                cls = "ICE-like"
            elif integ and hit["relaxase"]:
                cls = "IME-like"
            elif integ:
                cls = "integrative element"
            elif hit["relaxase"] or hit["MPF"]:
                cls = "MOB/MPF module"
            elif hit["transposase"] or tn_end:
                cls = "transposon acquisition"
            else:
                cls = "other"
            classes.append(cls)
            uniq = lambda lst: ";".join(dict.fromkeys(lst))
            seg_rows.append({
                "plasmid": x, "cluster": cl["cluster"], "length": xl, "start": a, "end": e_pos, "segment_bp": seg_bp,
                "class": cls, "replication": uniq(hit["replication"]), "relaxase": uniq(hit["relaxase"]),
                "MPF": uniq(hit["MPF"]), "integrase": uniq(hit["integrase"] + int_end),
                "transposase": uniq(hit["transposase"] + tn_end),
                "phage": uniq([lab for _, _, lab in pro] + hit["phage gene"]),
                "genes": sum(inside(g[0], g[1], region) for g in genes.get(x, [])),
                "best_match": donor or "", "best_match_cluster": clusters.get(donor, {}).get("cluster", ""),
                "best_match_type": clusters.get(donor, {}).get("cluster_type", ""),
                "best_match_fraction": round(dcov / seg_bp, 3), "source_in_dataset": int(in_data)})
            gene_rows += [{"type": f"segment end ({cls})", "plasmid": x, "start": g[0], "end": g[1], "locus": g[2],
                           "gene": g[3], "product": g[4]} for g in ends if MOBILE.search(g[4])]

        summary_rows.append({
            "plasmid": x, "cluster": cl["cluster"], "length": xl, "parents": len(parents),
            "parent_backbone_clusters": "|".join(systems[s]["backbone_cluster"] for s, _, _ in parents),
            "covered_by_parents_bp": length(core), "extra_segments": len(classes),
            "extra_bp": xl - length(core), "cointegrate": int(bool(best_pair)),
            **{c: classes.count(c) for c in CLASSES}})

    def write(path, rows, cols):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    write("compound_parents.tsv", parent_rows,
          ["plasmid", "cluster", "length", "circular", "systems", "plasmid_parents", "non_plasmid_backbones", "note"])
    write("extra_segments.tsv", seg_rows,
          ["plasmid", "cluster", "length", "start", "end", "segment_bp", "class", "replication", "relaxase", "MPF",
           "integrase", "transposase", "phage", "genes", "best_match", "best_match_cluster", "best_match_type",
           "best_match_fraction", "source_in_dataset"])
    write("cointegrates.tsv", coint_rows,
          ["plasmid", "cluster", "length", "system_1", "backbone_cluster_1", "backbone_1", "covered_kb_1", "fraction_1",
           "system_2", "backbone_cluster_2", "backbone_2", "covered_kb_2", "fraction_2", "overlap", "union",
           "junctions", "mobile_genes_at_junctions"])
    write("integration_summary.tsv", summary_rows,
          ["plasmid", "cluster", "length", "parents", "parent_backbone_clusters", "covered_by_parents_bp",
           "extra_segments", "extra_bp", "cointegrate"] + CLASSES)
    write("cointegrate_junction_genes.tsv", gene_rows, ["type", "plasmid", "start", "end", "locus", "gene", "product"])
    plot(summary_rows, seg_rows)
    n_par = sum(1 for r in parent_rows if r["plasmid_parents"])
    print(f"{len(parent_rows)} compound accessory plasmids, {n_par} with a plasmid parent; {len(seg_rows)} extra "
          f"segments >= {args.min_segment} bp; cointegrates: {len(coint_rows)}")


def plot(summary_rows, seg_rows):
    """One bar per compound plasmid (by length): parent-covered sequence and extra segments by class."""
    os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    rows = sorted(summary_rows, key=lambda r: r["length"])
    fig, ax = plt.subplots(figsize=(7, max(3, 0.06 * len(rows) + 1)))
    segs = defaultdict(list)
    for s in seg_rows:
        segs[s["plasmid"]].append(s)
    for i, r in enumerate(rows):
        ax.barh(i, r["length"] / 1000, color=COLORS["parent"], height=0.8, linewidth=0)
        for s in segs[r["plasmid"]]:
            for a, e in unwrap(s["start"], s["start"] + s["segment_bp"] - 1, r["length"]):
                ax.barh(i, (e - a + 1) / 1000, left=(a - 1) / 1000, color=COLORS[s["class"]], height=0.8, linewidth=0)
        if r["cointegrate"]:
            ax.text(r["length"] / 1000 + 3, i, "cointegrate", va="center", fontsize=5)
    used = {s["class"] for s in seg_rows}
    ax.set_yticks([])
    ax.set_ylabel(f"Compound accessory plasmids with a plasmid parent (n = {len(rows)})")
    ax.set_xlabel("Position (kb)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=[Patch(color=COLORS["parent"], label="covered by parent plasmid(s)")] +
              [Patch(color=COLORS[c], label=c) for c in CLASSES if c in used],
              frameon=False, fontsize=6, loc="lower right")
    fig.savefig("integration_map.png", dpi=300, bbox_inches="tight")
    fig.savefig("integration_map.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
