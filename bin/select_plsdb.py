#!/usr/bin/env python3
"""
Select PLSDB plasmids related to the target species.

A PLSDB record is selected when:
  1. metadata: the keyword (e.g. 'meliloti') appears, case-insensitive, in the
     MOB-typer `mash_neighbor_identification` or in `NUCCORE_Description`; or
  2. mash (optional): its accession is a hit in any `mash dist` table given.

PLSDB nuccore.csv (record metadata) and typing.csv (MOB-typer results) are
joined on NUCCORE_ACC; every nuccore record must have exactly one typing row.

Outputs:
  <prefix>.ids.txt       selected accessions, one per line (for seqkit grep)
  <prefix>.sources.tsv   accession <tab> source(s) (metadata, mash or metadata;mash)
  <prefix>.circular.tsv  accession <tab> 1 (circular) / 0 (linear), MobMess --complete input
"""
import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nuccore", required=True, help="PLSDB nuccore.csv")
    p.add_argument("--typing", required=True, help="PLSDB typing.csv (MOB-typer)")
    p.add_argument("--keyword", required=True, help="case-insensitive keyword, e.g. meliloti")
    p.add_argument("--mash", nargs="*", default=[], help="mash dist tables (reference ID in column 1)")
    p.add_argument("--prefix", default="plsdb_related")
    return p.parse_args()


def check_columns(reader, required, path):
    missing = required - set(reader.fieldnames or [])
    if missing:
        sys.exit(f"ERROR: {path} is missing columns: {', '.join(sorted(missing))}")


def main():
    args = parse_args()
    keyword = args.keyword.lower()
    csv.field_size_limit(sys.maxsize)

    # MOB-typer nearest-neighbour species per accession (typing.csv rows are not
    # in nuccore.csv order, so join on the accession, never by row)
    neighbour = {}
    with open(args.typing, newline="") as fh:
        reader = csv.DictReader(fh)
        check_columns(reader, {"NUCCORE_ACC", "mash_neighbor_identification"}, args.typing)
        for row in reader:
            acc = row["NUCCORE_ACC"]
            if acc in neighbour:
                sys.exit(f"ERROR: duplicated accession {acc} in {args.typing}")
            neighbour[acc] = row["mash_neighbor_identification"] or ""

    # Topology for every record, and keyword-based selection
    topology = {}
    sources = {}
    with open(args.nuccore, newline="") as fh:
        reader = csv.DictReader(fh)
        check_columns(reader, {"NUCCORE_ACC", "NUCCORE_Description", "NUCCORE_Topology"}, args.nuccore)
        for row in reader:
            acc = row["NUCCORE_ACC"]
            if acc in topology:
                sys.exit(f"ERROR: duplicated accession {acc} in {args.nuccore}")
            if acc not in neighbour:
                sys.exit(f"ERROR: {acc} from {args.nuccore} has no row in {args.typing} (different PLSDB versions?)")
            topology[acc] = row["NUCCORE_Topology"]
            text = f"{neighbour[acc]} {row['NUCCORE_Description'] or ''}".lower()
            if keyword in text:
                sources.setdefault(acc, set()).add("metadata")

    # Add every PLSDB accession hit by mash (thresholds already applied by mash dist -d/-v)
    for mash_file in args.mash:
        with open(mash_file) as fh:
            for line in fh:
                if line.strip():
                    acc = line.split("\t", 1)[0]
                    if acc not in topology:
                        sys.exit(f"ERROR: mash hit {acc} is not in {args.nuccore} (sketch/PLSDB version mismatch?)")
                    sources.setdefault(acc, set()).add("mash")

    if not sources:
        sys.exit(f"ERROR: no PLSDB record selected with keyword '{args.keyword}'")

    accs = sorted(sources)
    with open(f"{args.prefix}.ids.txt", "w") as ids, \
         open(f"{args.prefix}.sources.tsv", "w") as src, \
         open(f"{args.prefix}.circular.tsv", "w") as circ:
        src.write("accession\tsource\n")
        for acc in accs:
            ids.write(f"{acc}\n")
            src.write(f"{acc}\t{';'.join(sorted(sources[acc]))}\n")
            circ.write(f"{acc}\t{1 if topology[acc] == 'circular' else 0}\n")

    n_meta = sum("metadata" in s for s in sources.values())
    n_mash = sum("mash" in s for s in sources.values())
    print(f"Selected {len(accs)} PLSDB records (metadata: {n_meta}, mash: {n_mash})")


if __name__ == "__main__":
    main()
