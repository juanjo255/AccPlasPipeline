#!/usr/bin/env python3
"""
Keep CAZymes supported by at least --min-tools dbCAN tools (run_dbcan overview.tsv).

The publication required more than one tool (`#ofTools > 1`). Its DIAMOND search
returned nothing, so in practice dbCAN-HMM and dbCAN-sub both had to agree; run
dbCAN with `--methods hmm,dbCANsub` to reproduce that exactly.

Outputs:
  cazymes.tsv  plasmid, locus_tag, family, recommended result, substrate, tools
  cazymes.faa  their protein sequences (input of the optional SignalP step)
"""
import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--overview", required=True, help="run_dbcan overview.tsv")
    p.add_argument("--proteins", required=True, help="proteins given to dbCAN (headers <locus>__<plasmid>)")
    p.add_argument("--min-tools", type=int, default=2)
    return p.parse_args()


def main():
    args = parse_args()
    kept = []
    with open(args.overview, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["#ofTools"]) >= args.min_tools:
                locus, _, plasmid = r["Gene ID"].partition("__")
                rec = r["Recommend Results"]
                kept.append({"plasmid": plasmid, "locus_tag": locus, "protein_id": r["Gene ID"],
                             "family": rec.split("_")[0] if rec not in ("", "-") else "-",
                             "recommended": rec, "substrate": r.get("Substrate", "-"),
                             "dbCAN_hmm": r["dbCAN_hmm"], "dbCAN_sub": r["dbCAN_sub"],
                             "DIAMOND": r["DIAMOND"], "n_tools": r["#ofTools"]})

    cols = ["plasmid", "locus_tag", "protein_id", "family", "recommended", "substrate",
            "dbCAN_hmm", "dbCAN_sub", "DIAMOND", "n_tools"]
    with open("cazymes.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept)

    ids = {k["protein_id"] for k in kept}
    found, write = 0, False
    with open(args.proteins) as fh, open("cazymes.faa", "w") as out:
        for line in fh:
            if line.startswith(">"):
                write = line[1:].split(None, 1)[0] in ids
                found += write
            if write:
                out.write(line)
    if found != len(ids):
        sys.exit(f"ERROR: {len(ids) - found} CAZyme proteins not found in {args.proteins}")
    print(f"{len(kept)} CAZymes (>= {args.min_tools} tools) in {len({k['plasmid'] for k in kept})} plasmids")


if __name__ == "__main__":
    main()
