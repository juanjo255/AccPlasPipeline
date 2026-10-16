#!/usr/bin/env python3
"""
RepABC operons and one representative RepA per plasmid (publication rules).

From Bakta features whose product matches REP_REGEX (RepA/RepB/RepC), for each
repA (in order of sequence and start) the window of 4 kb downstream (+ strand)
or upstream (- strand) is collected. If the window holds repA, repB and repC,
that operon is complete: it replaces whatever was collected for the sequence and
later repA of the sequence are skipped. Otherwise the window is added.
A sequence gets a representative RepA only when exactly one repA row remains
(a complete operon, or a single repA).

With --wrap-circular, windows of circular sequences also wrap around the origin.
The publication notebook had this code but it never ran (circularity was
compared as text), so it is off by default.

Outputs:
  <prefix>.operons.tsv  rep features kept per sequence (with `complete_operon`)
  <prefix>.summary.tsv  sequence, rep genes combination, complete operon, RepA locus
  <prefix>.repA.tsv     sequence <tab> representative RepA locus tag
"""
import argparse
import csv
import re
import sys

# Same pattern as the publication notebook (the comma inside the class is literal)
REP_REGEX = re.compile(r"[rR]ep[Aa,Bb,Cc]")
REP_PARTITION_REGEX = re.compile(
    r"(?:^(?:(?:plasmid[\s/_-]*)?replication)|partition|para\s|parb\s)", re.IGNORECASE
)
BAKTA_COLUMNS = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]
WINDOW = 4000


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bakta", required=True, help="Bakta <prefix>.tsv")
    p.add_argument("--clusters", required=True, help="table with contig, circular, length (MobMess filtered clusters)")
    p.add_argument("--ids", help="only these sequences (one per line), e.g. plasmids kept in Step 6")
    p.add_argument("--wrap-circular", action="store_true", help="let windows wrap around circular sequences")
    p.add_argument("--prefix", default="repabc")
    return p.parse_args()


def main():
    args = parse_args()
    keep = None
    if args.ids:
        with open(args.ids) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    with open(args.clusters, newline="") as fh:
        info = {r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")}

    # Rep/partition features first (as in the notebook), then RepA/B/C among them
    reps = []
    with open(args.bakta) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            row = dict(zip(BAKTA_COLUMNS, line.rstrip("\n").split("\t")))
            if keep is not None and row["Sequence Id"] not in keep:
                continue
            product = row.get("Product", "")
            if REP_PARTITION_REGEX.search(product) and REP_REGEX.search(product):
                row["Start"], row["Stop"] = int(row["Start"]), int(row["Stop"])
                reps.append(row)
    reps.sort(key=lambda r: (r["Sequence Id"], r["Start"]))

    by_seq = {}
    for row in reps:
        by_seq.setdefault(row["Sequence Id"], []).append(row)

    complete, collected, n_unnamed = set(), {}, 0
    for row in reps:
        if row["Gene"] != "repA":
            continue
        seq = row["Sequence Id"]
        if seq in complete:
            continue
        rows = by_seq[seq]
        start, stop = row["Start"], row["Stop"]
        circular = info.get(seq, {}).get("circular") == "1"
        length = int(info[seq]["length"]) if seq in info else None
        if row["Strand"] == "+":
            window = [r for r in rows if r["Start"] >= start and r["Stop"] < start + WINDOW]
            if args.wrap_circular and circular and length and start + WINDOW > length:
                window += [r for r in rows if r["Stop"] < start + WINDOW - length]
        else:
            window = [r for r in rows if r["Start"] >= start - WINDOW and r["Stop"] <= stop]
            if args.wrap_circular and circular and length and start - WINDOW < 0:
                window += [r for r in rows if r["Start"] > length - (WINDOW - start)]

        # Deduplicate by locus tag within this window only
        seen, dedup = set(), []
        for r in window:
            if r["Locus Tag"] not in seen:
                seen.add(r["Locus Tag"])
                dedup.append(dict(r))
        genes = {r["Gene"] for r in dedup if r["Gene"]}
        n_unnamed += sum(1 for r in dedup if not r["Gene"])
        is_complete = "_".join(sorted(genes)) == "repA_repB_repC"
        for r in dedup:
            r["complete_operon"] = int(is_complete)
        if is_complete:
            complete.add(seq)
            collected[seq] = dedup
        else:
            collected.setdefault(seq, []).extend(dedup)

    if n_unnamed:
        print(f"WARNING: {n_unnamed} rep features without a gene name were ignored when naming operons", file=sys.stderr)

    cols = BAKTA_COLUMNS + ["complete_operon"]
    rep_a = {}
    with open(f"{args.prefix}.operons.tsv", "w", newline="") as out, \
         open(f"{args.prefix}.summary.tsv", "w", newline="") as summ:
        ops = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        ops.writeheader()
        sw = csv.writer(summ, delimiter="\t", lineterminator="\n")
        sw.writerow(["sequence", "rep_genes", "complete_operon", "repA_locus"])
        for seq, rows in collected.items():
            ops.writerows(rows)
            repa_rows = [r for r in rows if r["Gene"] == "repA"]
            locus = repa_rows[0]["Locus Tag"] if len(repa_rows) == 1 else "-"
            if locus != "-":
                rep_a[seq] = locus
            genes = "_".join(sorted({r["Gene"] for r in rows if r["Gene"]}))
            sw.writerow([seq, genes, int(seq in complete), locus])

    with open(f"{args.prefix}.repA.tsv", "w") as out:
        for seq, locus in rep_a.items():
            out.write(f"{seq}\t{locus}\n")

    print(f"{len(collected)} sequences with RepA, {len(complete)} with a complete RepABC operon, "
          f"{len(rep_a)} with one representative RepA")


if __name__ == "__main__":
    main()
