#!/usr/bin/env python3
"""
Circularity and read depth of the length-filtered contigs from Flye's assembly_info.txt.

Contig IDs are renamed exactly as in the contig length filter: any existing
`_<sample>` suffix is stripped, then re-added if --add-suffix is true. This
handles both raw Flye output and assembly_info files already renamed.

Outputs (only for contigs present in the filtered FASTA):
  <sample>.circular.tsv  contig <tab> 1 (circ. = Y) / 0 (circ. = N)
  <sample>.coverage.tsv  contig, depth (cov.), chromosome_depth, rel_coverage
                         (chromosome = largest contig of the whole assembly)
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

    circular, depth = {}, {}
    chrom_len, chrom_depth = -1, None
    with open(args.info, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"#seq_name", "length", "cov.", "circ."}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            sys.exit(f"ERROR: {args.info} does not look like a Flye assembly_info.txt (needs columns {sorted(required)})")
        for row in reader:
            if int(row["length"]) > chrom_len:
                chrom_len, chrom_depth = int(row["length"]), int(row["cov."])
            depth[rename(row["#seq_name"])] = int(row["cov."])
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

    with open(f"{args.sample}.coverage.tsv", "w") as out:
        out.write("contig\tdepth\tchromosome_depth\trel_coverage\n")
        for contig in wanted:
            rel = depth[contig] / chrom_depth if chrom_depth else float("nan")
            out.write(f"{contig}\t{depth[contig]}\t{chrom_depth}\t{rel:.3f}\n")

    n_circ = sum(circular[c] == "1" for c in wanted)
    print(f"[{args.sample}] {len(wanted)} contigs, {n_circ} circular")


if __name__ == "__main__":
    main()
