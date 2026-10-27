#!/usr/bin/env python3
"""
BGC novelty of the accessory plasmids (after Pulliam et al. 2026, bioRxiv 10.64898/2026.08.08.743306).

Per BGC (antiSMASH region):
  label             known = its BiG-SCAPE mix-bin GCF at --bigscape-cutoff (default 0.3) contains a
                    MIBiG BGC; novel otherwise
  antismashdb_match distance to the nearest antiSMASH-DB GCF (bigslice --query) <= --match-dist
                    (previously observed in public genomes, not characterized chemistry)
Per BiG-SLiCE GCF (own clustering, threshold 0.65; BGCs beyond it are singletons):
  novelty_mean_asdb4_dist  mean member distance to antiSMASH-DB GCFs (the paper's novelty score)
  antismashdb_match        mean <= --match-dist; frac_members_antismashdb_match; size_ge<N>

Outputs:
  bgc_novelty.tsv              one row per BGC
  gcf_novelty_bigslice.tsv     one row per BiG-SLiCE GCF or singleton (size >= N first, most novel first)
"""
import argparse
import csv
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from statistics import mean, median

MIBIG_RE = re.compile(r"^BGC\d{7}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bigscape-db", required=True, help="BiG-SCAPE 2 SQLite database")
    p.add_argument("--bigslice-db", required=True, help="BiG-SLiCE clustering result/data.db")
    p.add_argument("--query-db", required=True, help="BiG-SLiCE query report data.db (reports/<id>/data.db)")
    p.add_argument("--asdb-db", required=True, help="antiSMASH-DB BiG-SLiCE result/data.db (reference run)")
    p.add_argument("--bigscape-cutoff", type=float, default=0.3)
    p.add_argument("--match-dist", type=float, default=0.65)
    p.add_argument("--min-gcf-size", type=int, default=5)
    return p.parse_args()


def key_from_path(path):
    # Join key shared by both tools: GBK filename without extension
    return re.sub(r"\.gbk$", "", os.path.basename(path))


def connect_ro(path):
    return sqlite3.connect(f"file:{os.path.abspath(path)}?mode=ro", uri=True)


def load_bigscape(db_path, cutoff):
    con = connect_ro(db_path)
    run_id = con.execute("SELECT max(id) FROM run").fetchone()[0]
    rec = {}  # record_id -> (key, product, is_mibig)
    # Only region-level records get families (the DB also holds protocluster/cand_cluster/proto_core)
    for rid, path, product in con.execute(
            "SELECT r.id, g.path, r.product FROM bgc_record r JOIN gbk g ON g.id = r.gbk_id "
            "WHERE r.record_type = 'region'"):
        k = key_from_path(path)
        rec[rid] = (k, product, bool(MIBIG_RE.match(k)))
    n_mibig = sum(v[2] for v in rec.values())
    if n_mibig == 0:
        # Without MIBiG records every BGC would silently be "novel"
        sys.exit(f"ERROR: no MIBiG records in {db_path}; BiG-SCAPE must run with the MIBiG GenBank files")

    fam = {}
    for rid, fid, c in con.execute(
            "SELECT rf.record_id, f.id, f.cutoff FROM bgc_record_family rf JOIN family f ON f.id = rf.family_id "
            "WHERE f.bin_label = 'mix' AND f.run_id = ?", (run_id,)):
        if rid in rec and abs(c - cutoff) < 1e-9:
            fam[rid] = fid
    if not fam:
        sys.exit(f"ERROR: no mix-bin families at cutoff {cutoff} in {db_path} (run BiG-SCAPE with --mix)")
    fam_mibig = defaultdict(set)
    for rid, fid in fam.items():
        if rec[rid][2]:
            fam_mibig[fid].add(rec[rid][0])

    # Closest MIBiG BGC by BiG-SCAPE distance (pairs never compared count as 1.0)
    mibig_ids = {rid for rid, v in rec.items() if v[2]}
    best = {}
    for a, b, d in con.execute("SELECT record_a_id, record_b_id, distance FROM distance"):
        for q, r in ((a, b), (b, a)):
            if r in mibig_ids and q in rec and q not in mibig_ids and d < best.get(q, (1.0, ""))[0]:
                best[q] = (d, rec[r][0])
    con.close()

    out = {}
    for rid, (k, product, is_mibig) in rec.items():
        if is_mibig:
            continue
        fid = fam.get(rid)
        d, near = best.get(rid, (1.0, ""))
        out[k] = {"bigscape_product": product, "bigscape_gcf": fid if fid is not None else "",
                  "mibig_in_gcf": ";".join(sorted(fam_mibig.get(fid, ()))),
                  "nearest_mibig": near, "nearest_mibig_dist": round(d, 4)}
    print(f"BiG-SCAPE: {len(out)} query regions, {n_mibig} MIBiG regions", file=sys.stderr)
    return out


def load_bigslice_own(db_path):
    con = connect_ro(db_path)
    thr = con.execute("SELECT threshold FROM clustering ORDER BY id DESC LIMIT 1").fetchone()[0]
    out = {}
    for fname, edge, gcf, val in con.execute(
            "SELECT b.orig_filename, b.on_contig_edge, g.id_in_run, m.membership_value FROM gcf_membership m "
            "JOIN bgc b ON b.id = m.bgc_id JOIN gcf g ON g.id = m.gcf_id "
            "WHERE m.rank = 0 AND g.clustering_id = (SELECT max(id) FROM clustering)"):
        k = key_from_path(fname)
        # Beyond the clustering threshold a BGC is not a member: its own group (nearest GCF kept)
        out[k] = {"bigslice_gcf": str(gcf) if val <= thr else f"singleton:{k}",
                  "bigslice_nearest_gcf": gcf, "bigslice_gcf_dist": round(val, 4),
                  "on_contig_edge": bool(edge)}
    con.close()
    return out, thr


def load_bigslice_query(query_db, asdb_db):
    src = connect_ro(asdb_db)
    thr = src.execute("SELECT threshold FROM clustering ORDER BY id DESC LIMIT 1").fetchone()[0]
    gcf_label = dict(src.execute("SELECT id, id_in_run FROM gcf").fetchall())
    src.close()
    out = {}
    q = connect_ro(query_db)
    for fname, gcf_id, val, rank in q.execute(
            "SELECT b.orig_filename, m.gcf_id, m.membership_value, m.rank FROM gcf_membership m "
            "JOIN bgc b ON b.id = m.bgc_id ORDER BY m.rank"):
        d = out.setdefault(key_from_path(fname), {"asdb4_hits": []})
        d["asdb4_hits"].append(f"{gcf_label.get(gcf_id, gcf_id)}:{val:.4f}")
        if rank == 0:
            d["asdb4_best_gcf"] = gcf_label.get(gcf_id, gcf_id)
            d["asdb4_best_dist"] = round(val, 4)
    q.close()
    for d in out.values():
        d["asdb4_hits"] = ";".join(d["asdb4_hits"])
    return out, thr


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    a = parse_args()
    bs = load_bigscape(a.bigscape_db, a.bigscape_cutoff)
    own, own_thr = load_bigslice_own(a.bigslice_db)
    qry, asdb_thr = load_bigslice_query(a.query_db, a.asdb_db)

    rows = []
    for k in sorted(set(bs) | set(own) | set(qry)):
        r = {"bgc": k, **bs.get(k, {}), **own.get(k, {}), **qry.get(k, {})}
        r["label"] = "known" if r.get("mibig_in_gcf") else "novel"
        d = r.get("asdb4_best_dist")
        r["antismashdb_match"] = "" if d is None else d <= a.match_dist
        rows.append(r)

    by_gcf = defaultdict(list)
    for r in rows:
        if "bigslice_gcf" in r:
            by_gcf[r["bigslice_gcf"]].append(r)
    size_col = f"size_ge{a.min_gcf_size}"
    gcf_rows = []
    for g, members in by_gcf.items():
        d = sorted(m["asdb4_best_dist"] for m in members if m.get("asdb4_best_dist") is not None)
        n_known = sum(m["label"] == "known" for m in members)
        gcf_rows.append({
            "bigslice_gcf": g, "n_bgcs": len(members), size_col: len(members) >= a.min_gcf_size,
            "novelty_mean_asdb4_dist": round(mean(d), 4) if d else "",
            "antismashdb_match": (mean(d) <= a.match_dist) if d else "",
            "frac_members_antismashdb_match": round(sum(m["antismashdb_match"] is True for m in members) / len(members), 3),
            "min_asdb4_dist": d[0] if d else "", "median_asdb4_dist": round(median(d), 4) if d else "",
            "max_asdb4_dist": d[-1] if d else "",
            "n_known": n_known, "n_novel": len(members) - n_known,
            "n_contig_edge": sum(bool(m.get("on_contig_edge")) for m in members),
            "products": ";".join(sorted({str(m.get("bigscape_product", "")) for m in members})),
            "mibig": ";".join(sorted({x for m in members for x in str(m.get("mibig_in_gcf", "")).split(";") if x})),
            "members": ";".join(m["bgc"] for m in members),
        })
    # GCFs meeting the size criterion first, then most novel first
    gcf_rows.sort(key=lambda x: (not x[size_col], x["novelty_mean_asdb4_dist"] == "",
                                 -(x["novelty_mean_asdb4_dist"] or 0)))

    bgc_cols = ["bgc", "bigscape_product", "label", "bigscape_gcf", "mibig_in_gcf", "nearest_mibig",
                "nearest_mibig_dist", "bigslice_gcf", "bigslice_nearest_gcf", "bigslice_gcf_dist",
                "on_contig_edge", "asdb4_best_gcf", "asdb4_best_dist", "antismashdb_match", "asdb4_hits"]
    write_tsv("bgc_novelty.tsv", rows, bgc_cols)
    write_tsv("gcf_novelty_bigslice.tsv", gcf_rows, list(gcf_rows[0]) if gcf_rows else ["bigslice_gcf"])

    counts = Counter(r["label"] for r in rows)
    print(f"BiG-SCAPE cutoff {a.bigscape_cutoff}; BiG-SLiCE threshold {own_thr}; antiSMASH-DB threshold {asdb_thr}",
          file=sys.stderr)
    print(f"BGCs {len(rows)}: known {counts['known']}, novel {counts['novel']}, "
          f"antismashdb_match {sum(r['antismashdb_match'] is True for r in rows)}; "
          f"GCF groups {len(gcf_rows)} ({sum(g[size_col] for g in gcf_rows)} with >= {a.min_gcf_size} BGCs)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
