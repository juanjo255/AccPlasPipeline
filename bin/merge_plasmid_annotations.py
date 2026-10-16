#!/usr/bin/env python3
"""
Merge per-plasmid Bakta outputs into one set of files for the annotation tools.

Outputs:
  plasmids.tsv           all Bakta features (one header block)
  plasmids_proteins.faa  all proteins, headers `<locus>__<plasmid> <product>`
                         (TASmania and dbCAN keep the plasmid in the protein ID)
  plasmids.gbff          all GenBank records (antiSMASH input)
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsv", nargs="+", required=True, help="per-plasmid Bakta .tsv")
    p.add_argument("--faa", nargs="+", required=True, help="per-plasmid Bakta .faa")
    p.add_argument("--gbff", nargs="+", required=True, help="per-plasmid Bakta .gbff")
    return p.parse_args()


def main():
    args = parse_args()

    # Features: keep Bakta's comment header once, sort files for a stable order
    locus_to_plasmid = {}
    with open("plasmids.tsv", "w") as out:
        for i, path in enumerate(sorted(args.tsv)):
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        if i == 0:
                            out.write(line)
                        continue
                    out.write(line)
                    f = line.rstrip("\n").split("\t")
                    if len(f) > 5 and f[5]:
                        if f[5] in locus_to_plasmid and locus_to_plasmid[f[5]] != f[0]:
                            sys.exit(f"ERROR: locus tag {f[5]} found in {locus_to_plasmid[f[5]]} and {f[0]}")
                        locus_to_plasmid[f[5]] = f[0]

    n_prot = 0
    with open("plasmids_proteins.faa", "w") as out:
        for path in sorted(args.faa):
            with open(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        locus, _, rest = line[1:].rstrip("\n").partition(" ")
                        if locus not in locus_to_plasmid:
                            sys.exit(f"ERROR: protein {locus} of {path} not in the Bakta TSVs")
                        out.write(f">{locus}__{locus_to_plasmid[locus]} {rest}\n")
                        n_prot += 1
                    else:
                        out.write(line)

    with open("plasmids.gbff", "w") as out:
        for path in sorted(args.gbff):
            with open(path) as fh:
                out.write(fh.read())

    print(f"{len(args.tsv)} plasmids, {n_prot} proteins merged")


if __name__ == "__main__":
    main()
