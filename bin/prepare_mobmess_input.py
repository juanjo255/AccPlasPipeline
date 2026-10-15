#!/usr/bin/env python3
"""
Pool sample contigs and PLSDB plasmids into one MobMess input.

Sequences are sorted by ID so the MobMess input (and its cluster numbering)
does not depend on the order in which samples finish. Fails on duplicated IDs
and when FASTA IDs and circularity rows do not match one-to-one.

Outputs: mobmess_input.fasta, mobmess_input.circular.tsv
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", nargs="+", required=True, help="FASTA files to pool")
    p.add_argument("--circular", nargs="+", required=True, help="contig <tab> 1/0 tables to pool")
    return p.parse_args()


def read_fasta(path, seqs):
    # Intent: collect header line and sequence lines per ID, rejecting duplicates
    seq_id = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].split(None, 1)[0]
                if seq_id in seqs:
                    sys.exit(f"ERROR: duplicated sequence ID {seq_id} (in {path})")
                seqs[seq_id] = [line if line.endswith("\n") else line + "\n"]
            elif seq_id is not None and line.strip():
                seqs[seq_id].append(line if line.endswith("\n") else line + "\n")


def main():
    args = parse_args()

    seqs = {}
    for path in args.fasta:
        read_fasta(path, seqs)

    circular = {}
    for path in args.circular:
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                contig, flag = line.rstrip("\n").split("\t")[:2]
                if contig in circular:
                    sys.exit(f"ERROR: duplicated circularity row for {contig} (in {path})")
                if flag not in ("0", "1"):
                    sys.exit(f"ERROR: circularity for {contig} must be 0 or 1, got '{flag}'")
                circular[contig] = flag

    no_flag = sorted(set(seqs) - set(circular))
    no_seq = sorted(set(circular) - set(seqs))
    if no_flag or no_seq:
        sys.exit("ERROR: FASTA and circularity tables do not match. "
                 f"Without circularity: {no_flag[:5]}; without sequence: {no_seq[:5]}")

    with open("mobmess_input.fasta", "w") as fa, open("mobmess_input.circular.tsv", "w") as tsv:
        for seq_id in sorted(seqs):
            fa.writelines(seqs[seq_id])
            tsv.write(f"{seq_id}\t{circular[seq_id]}\n")

    n_circ = sum(v == "1" for v in circular.values())
    print(f"MobMess input: {len(seqs)} sequences, {n_circ} circular")


if __name__ == "__main__":
    main()
