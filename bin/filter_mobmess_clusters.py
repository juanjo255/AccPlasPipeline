#!/usr/bin/env python3
"""
Keep MobMess clusters that look like real plasmids (rules of the publication).

Per cluster type:
  fragment               discarded
  backbone               all members < max_len
  compound               all members < max_len and >= 1 circular member
  maximal_not_in_system  all members < max_len, >= 1 circular member and
                         >= 1 member that is not a PLSDB reference

PLSDB-only backbone and compound clusters are kept on purpose: MobMess places
them in plasmid systems, so their relationship to the target plasmids is
established. A PLSDB-only maximal cluster is related only through the PLSDB
annotation used to select it, so it is not kept.

Outputs:
  <prefix>.filtered_clusters.tsv  MobMess clusters table rows that passed, plus
                                  `length` and `source` (sample / plsdb)
  <prefix>.filtered_ids.txt       contig IDs that passed
  <prefix>.filtered.fasta         sequences that passed
"""
import argparse
import collections
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clusters", required=True, help="MobMess <prefix>_clusters.txt")
    p.add_argument("--fasta", required=True, help="FASTA given to MobMess (for lengths)")
    p.add_argument("--plsdb-ids", required=True, help="PLSDB accessions in the MobMess input")
    p.add_argument("--max-len", type=int, default=1000000, help="members must be shorter than this (bp)")
    p.add_argument("--prefix", default="mobmess")
    return p.parse_args()


def fasta_lengths(path):
    # Intent: sequence length per ID (first word of the header)
    lengths, seq_id, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    lengths[seq_id] = n
                seq_id, n = line[1:].split(None, 1)[0], 0
            else:
                n += len(line.strip())
    if seq_id is not None:
        lengths[seq_id] = n
    return lengths


def write_fasta(path, keep, out_path):
    # Intent: stream the input FASTA, writing only records whose ID is kept
    with open(path) as fh, open(out_path, "w") as out:
        write = False
        for line in fh:
            if line.startswith(">"):
                write = line[1:].split(None, 1)[0] in keep
            if write:
                out.write(line)


def keep_cluster(cluster_type, members, max_len):
    if cluster_type == "fragment":
        return False
    if not all(m["length"] < max_len for m in members):
        return False
    if cluster_type == "backbone":
        return True
    has_circular = any(m["circular"] == "1" for m in members)
    if cluster_type == "compound":
        return has_circular
    if cluster_type == "maximal_not_in_system":
        return has_circular and any(m["source"] == "sample" for m in members)
    sys.exit(f"ERROR: unknown MobMess cluster type '{cluster_type}'")


def main():
    args = parse_args()
    lengths = fasta_lengths(args.fasta)
    plsdb = {l.strip() for l in open(args.plsdb_ids) if l.strip()}

    with open(args.clusters, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames + ["length", "source"]
        rows = list(reader)

    clusters = collections.defaultdict(list)
    for row in rows:
        if row["contig"] not in lengths:
            sys.exit(f"ERROR: {row['contig']} is in {args.clusters} but not in {args.fasta}")
        row["length"] = lengths[row["contig"]]
        row["source"] = "plsdb" if row["contig"] in plsdb else "sample"
        clusters[row["cluster"]].append(row)

    kept = []
    for members in clusters.values():
        types = {m["cluster_type"] for m in members}
        if len(types) != 1:
            sys.exit(f"ERROR: cluster {members[0]['cluster']} has several types: {types}")
        if keep_cluster(types.pop(), members, args.max_len):
            kept.extend(members)

    with open(f"{args.prefix}.filtered_clusters.tsv", "w", newline="") as out, \
         open(f"{args.prefix}.filtered_ids.txt", "w") as ids:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in kept:
            writer.writerow(row)
            ids.write(f"{row['contig']}\n")

    write_fasta(args.fasta, {r["contig"] for r in kept}, f"{args.prefix}.filtered.fasta")

    n_clusters = len({r["cluster"] for r in kept})
    n_sample = sum(r["source"] == "sample" for r in kept)
    print(f"Kept {len(kept)} sequences ({n_sample} sample, {len(kept) - n_sample} PLSDB) in {n_clusters} clusters")


if __name__ == "__main__":
    main()
