#!/usr/bin/env python3
"""
Circularity of the length-filtered contigs from Flye's assembly_info.txt.

Contig IDs are renamed exactly as in the contig length filter: any existing
`_<sample>` suffix is stripped, then re-added if --add-suffix is true. This
handles both raw Flye output and assembly_info files already renamed.

Output: <sample>.circular.tsv  contig <tab> 1 (circ. = Y) / 0 (circ. = N),
        only for contigs present in the filtered FASTA.
"""
import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--info", required=True, help="Flye assembly_info.txt")
    p.add_argument("--fasta", required=True, help="length-filtered contigs FASTA")
    p.add_argument("--sample", required=True)
    p.add_argument("--add-suffix", choices=["true", "false"], required=True)
    return p.parse_args()


def main():
    args = parse_args()
    suffix = f"_{args.sample}"

    def rename(contig):
        base = contig[: -len(suffix)] if contig.endswith(suffix) else contig
        return base + suffix if args.add_suffix == "true" else base

    with open(args.fasta) as fh:
        wanted = [line[1:].split(None, 1)[0] for line in fh if line.startswith(">")]

    circular = {}
    with open(args.info, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "#seq_name" or "circ." not in reader.fieldnames:
            sys.exit(f"ERROR: {args.info} does not look like a Flye assembly_info.txt (needs '#seq_name' and 'circ.' columns)")
        for row in reader:
            flag = row["circ."]
            if flag not in ("Y", "N"):
                sys.exit(f"ERROR: unexpected circ. value '{flag}' for {row['#seq_name']} in {args.info}")
            circular[rename(row["#seq_name"])] = "1" if flag == "Y" else "0"

    missing = [c for c in wanted if c not in circular]
    if missing:
        sys.exit(f"ERROR: [{args.sample}] contigs not found in {args.info}: {', '.join(missing[:5])}")

    with open(f"{args.sample}.circular.tsv", "w") as out:
        for contig in wanted:
            out.write(f"{contig}\t{circular[contig]}\n")

    n_circ = sum(circular[c] == "1" for c in wanted)
    print(f"[{args.sample}] {len(wanted)} contigs, {n_circ} circular")


if __name__ == "__main__":
    main()
