#!/usr/bin/env python3
"""
Replicon groups from RepA clusters (Cavassim et al. 2020 naming).

`extract` writes the representative RepA proteins (one per plasmid, from
repabc_operons.py) as <locus>_<sequence> followed by the reference RepA
proteins, ready for `mmseqs easy-cluster`.

`name` names the MMseqs2 clusters: a cluster containing a reference RepA takes
its group name (e.g. Rh07); other clusters are named Si01, Si02, ... by
decreasing size (ties broken by representative ID). With --names (sequence <tab>
group, e.g. the publication's assignments), a cluster without a reference takes
the listed name held by most of its members; new clusters are numbered after
the highest Si listed. This does not depend on which member MMseqs2 picks as
representative.

Outputs:
  extract: <prefix>.faa, <prefix>.map.tsv (header <tab> sequence)
  name:    <prefix>.replicon_groups.tsv (sequence, repA_locus, replicon_group,
           cluster_representative, group_size)
"""
import argparse
import collections
import csv
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--repa", required=True, help="sequence <tab> RepA locus tag")
    e.add_argument("--faa", required=True, help="Bakta <prefix>.faa")
    e.add_argument("--references", required=True, help="reference RepA proteins (headers with group names, e.g. RepA_Rh01)")
    e.add_argument("--prefix", default="repA")
    n = sub.add_parser("name")
    n.add_argument("--clusters", required=True, help="mmseqs <prefix>_cluster.tsv (representative <tab> member)")
    n.add_argument("--map", required=True, help="<prefix>.map.tsv from extract")
    n.add_argument("--references", required=True, help="same reference FASTA as in extract")
    n.add_argument("--names", help="optional fixed names: sequence <tab> replicon group")
    n.add_argument("--prefix", default="repA")
    return p.parse_args()


def read_fasta(path):
    # Intent: {id: sequence} keeping file order
    seqs, seq_id = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].split(None, 1)[0]
                seqs[seq_id] = []
            elif seq_id is not None:
                seqs[seq_id].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def reference_groups(path):
    # Reference header -> group name (only reference headers may supply a group)
    groups = {}
    for ref in read_fasta(path):
        m = re.search(r"Rh[0-9]+", ref)
        if not m:
            sys.exit(f"ERROR: reference header '{ref}' has no group name (e.g. Rh01)")
        groups[ref] = m.group(0)
    return groups


def extract(args):
    proteins = read_fasta(args.faa)
    refs = read_fasta(args.references)
    reference_groups(args.references)
    with open(args.repa) as fh:
        repa = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    wanted = {locus: seq for seq, locus in repa}
    missing = [l for l in wanted if l not in proteins]
    if missing:
        sys.exit(f"ERROR: RepA loci not in {args.faa}: {', '.join(missing[:5])}")
    with open(f"{args.prefix}.faa", "w") as out, open(f"{args.prefix}.map.tsv", "w") as mp:
        # Plasmid RepA first (in .faa order), then references: input order sets MMseqs2 representatives
        for locus, prot in proteins.items():
            if locus in wanted:
                header = f"{locus}_{wanted[locus]}"
                out.write(f">{header}\n{prot}\n")
                mp.write(f"{header}\t{wanted[locus]}\t{locus}\n")
        for ref, prot in refs.items():
            out.write(f">{ref}\n{prot}\n")
    print(f"{len(wanted)} plasmid RepA + {len(refs)} reference RepA written")


def name(args):
    groups = reference_groups(args.references)
    with open(args.map) as fh:
        header_info = {h: (s, l) for h, s, l in (line.rstrip("\n").split("\t") for line in fh if line.strip())}
    members = collections.defaultdict(list)
    with open(args.clusters) as fh:
        for line in fh:
            if line.strip():
                rep, member = line.rstrip("\n").split("\t")
                members[rep].append(member)

    names = {}
    for rep, mem in members.items():
        rh = sorted({groups[m] for m in mem if m in groups})
        if len(rh) > 1:
            print(f"WARNING: cluster {rep} contains several reference groups {rh}; using {rh[0]}", file=sys.stderr)
        if rh:
            names[rep] = rh[0]
    fixed = {}
    if args.names:
        with open(args.names) as fh:
            fixed = dict(line.rstrip("\n").split("\t")[:2] for line in fh if line.strip() and not line.startswith("#"))
    # Majority of members' listed names (sequence IDs; Bakta locus tags can change between runs)
    for rep, mem in members.items():
        if rep in names:
            continue
        votes = collections.Counter(fixed[header_info[m][0]] for m in mem
                                    if m in header_info and header_info[m][0] in fixed)
        if votes:
            names[rep] = votes.most_common(1)[0][0]
    start = max([int(v[2:]) for v in fixed.values() if re.fullmatch(r"Si[0-9]+", v)] + [0])
    unnamed = sorted((r for r in members if r not in names), key=lambda r: (-len(members[r]), r))
    for i, rep in enumerate(unnamed, start + 1):
        names[rep] = f"Si{i:02d}"
    if len(set(names.values())) != len(names):
        dup = [n for n, c in collections.Counter(names.values()).items() if c > 1]
        sys.exit(f"ERROR: several clusters got the same replicon group name: {', '.join(dup)}")

    n_written = 0
    with open(f"{args.prefix}.replicon_groups.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["sequence", "repA_locus", "replicon_group", "cluster_representative", "group_size"])
        for rep, mem in sorted(members.items(), key=lambda kv: names[kv[0]]):
            for m in mem:
                if m in header_info:
                    seq, locus = header_info[m]
                    w.writerow([seq, locus, names[rep], rep, len(mem)])
                    n_written += 1
    print(f"{n_written} plasmids in {len({names[r] for r in members})} replicon groups "
          f"({sum(1 for r in names.values() if r.startswith('Rh'))} with a reference group)")


if __name__ == "__main__":
    a = parse_args()
    extract(a) if a.cmd == "extract" else name(a)
