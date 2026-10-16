#!/usr/bin/env python3
"""
Toxin-antitoxin (TA) genes like the TASer web server (TASmania HMMs).

TASer (https://shiny.bioinformatics.unibe.ch/apps/tasmania/) can only be run on its
web page, so its method is reproduced locally. The rules were inferred from the
publication's TASer reports (Prodigal training on the submitted sequences reproduced
100% of TASer's gene coordinates):
  - genes are predicted with Prodigal (pyrodigal, single mode) trained on all
    submitted plasmids together; gene IDs are <plasmid>_<n> as in TASer;
  - every protein is scored alone against the TASmania HMMs (E-values with Z = 1),
    so results do not depend on how many plasmids are analysed;
  - each gene keeps its best (lowest E-value) HMM and is reported at
    E < --report-evalue (0.01, HMMER's default inclusion threshold);
  - toxin (T) / antitoxin (A) comes from the HMM file name (T.<id>.hmm3 / A.<id>.hmm3);
  - pseudo-operons are runs of genes on the same strand with intergenic distance
    between -100 and +100 nt; the structure lists them in genomic order as T, A or
    x (no TA hit), e.g. "xTA".
Each TA gene is linked to the Bakta CDS it overlaps (>= 50% of both lengths, as in
the publication notebook) for its gene name and product.
Counted TA genes (`count` = 1) have E < --evalue (1e-3) and are not RepA/B/C proteins.

Outputs:
  <prefix>.genes.faa  predicted proteins (<plasmid>_<n>)
  <prefix>.hits.csv   every HMM hit
  <prefix>.ta.tsv     one row per reported TA gene
"""
import argparse
import csv
import glob
import os
import sys
import warnings

import pyhmmer
import pyrodigal
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile

REPABC = {"repA", "repB", "repC"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="plasmid nucleotide sequences")
    p.add_argument("--hmm-dir", required=True, help="directory with TASmania T.*.hmm3 / A.*.hmm3 files")
    p.add_argument("--bakta", required=True, help="Bakta TSV of the same plasmids (gene name, product)")
    p.add_argument("--z", type=float, default=1.0, help="database size for E-values (TASer: 1)")
    p.add_argument("--report-evalue", type=float, default=1e-2, help="report TA genes below this E-value (TASer-like)")
    p.add_argument("--evalue", type=float, default=1e-3, help="E-value cutoff for counting (publication: 1e-3)")
    p.add_argument("--min-gap", type=int, default=-100, help="min intergenic distance (overlap) within a pseudo-operon (nt)")
    p.add_argument("--max-gap", type=int, default=100, help="max intergenic distance within a pseudo-operon (nt)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--prefix", default="tasmania")
    return p.parse_args()


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                if name in seqs:
                    sys.exit(f"ERROR: duplicated sequence ID {name} in {path}")
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_hmms(hmm_dir):
    # T/A from the file name; the file stem is the TASmania HMM ID (e.g. T.17528)
    hmms, info = [], {}
    for path in sorted(glob.glob(os.path.join(hmm_dir, "*.hmm3"))):
        stem = os.path.basename(path)[: -len(".hmm3")]
        kind = stem.split(".", 1)[0]
        if kind not in ("T", "A"):
            sys.exit(f"ERROR: HMM file {path} is not named T.<id>.hmm3 or A.<id>.hmm3")
        with HMMFile(path) as fh:
            for hmm in fh:
                name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                info[name] = (stem, kind)
                hmms.append(hmm)
    if not hmms:
        sys.exit(f"ERROR: no *.hmm3 files in {hmm_dir}")
    return hmms, info


def read_bakta(path):
    cds = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and f[1] == "cds":
                cds.setdefault(f[0], []).append((int(f[2]), int(f[3]), f[5], f[6], f[7]))
    return cds


def bakta_match(cds, start, stop):
    # Bakta CDS overlapping >= 50% of both genes (publication notebook rule)
    hits = []
    for s, e, locus, gene, product in cds:
        ov = min(e, stop) - max(s, start) + 1
        if ov > 0 and ov / (stop - start + 1) >= 0.5 and ov / (e - s + 1) >= 0.5:
            hits.append((locus, gene, product))
    return hits


def main():
    args = parse_args()
    seqs = read_fasta(args.fasta)
    hmms, hmm_info = read_hmms(args.hmm_dir)
    bakta = read_bakta(args.bakta)

    # Prodigal single mode trained on all plasmids together (as TASer on a submission);
    # with too little sequence to train (Prodigal recommends >= 100 kb), use metagenome mode
    total = sum(len(s) for s in seqs.values())
    if total >= 100000:
        finder = pyrodigal.GeneFinder()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            finder.train(*[s.encode() for s in seqs.values()])
    else:
        print(f"WARNING: only {total} bp to train Prodigal; using metagenome mode (gene calls may differ from TASer)",
              file=sys.stderr)
        finder = pyrodigal.GeneFinder(meta=True)
    genes = {}
    alphabet = Alphabet.amino()
    proteins = []
    with open(f"{args.prefix}.genes.faa", "w") as faa:
        for plasmid, seq in seqs.items():
            for n, g in enumerate(finder.find_genes(seq.encode()), 1):
                gid = f"{plasmid}_{n}"
                strand = "+" if g.strand == 1 else "-"
                genes[gid] = (plasmid, n, g.begin, g.end, strand)
                prot = g.translate().rstrip("*")
                faa.write(f">{gid}\n{prot}\n")
                proteins.append(TextSequence(name=gid.encode(), sequence=prot).digitize(alphabet))
    targets = pyhmmer.easel.DigitalSequenceBlock(alphabet, proteins)

    rows = []
    for top_hits in pyhmmer.hmmsearch(hmms, targets, cpus=args.threads, Z=args.z):
        query = top_hits.query.name
        query = query.decode() if isinstance(query, bytes) else query
        hmm_id, kind = hmm_info[query]
        for hit in top_hits:
            name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            rows.append({"hmm_id": hmm_id, "ta_type": kind, "hmm_name": query, "hmm_length": top_hits.query.M,
                         "gene_id": name, "plasmid": genes[name][0], "bitscore": hit.score,
                         "evalue": hit.evalue, "length": hit.length})
    cols = ["hmm_id", "ta_type", "hmm_name", "hmm_length", "gene_id", "plasmid", "bitscore", "evalue", "length"]
    with open(f"{args.prefix}.hits.csv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # Best HMM per gene (ties: alphabetical HMM ID), reported below the TASer-like cutoff
    best = {}
    for r in rows:
        if r["evalue"] < args.report_evalue:
            g = r["gene_id"]
            if g not in best or (r["evalue"], r["hmm_id"]) < (best[g]["evalue"], best[g]["hmm_id"]):
                best[g] = r

    # Pseudo-operons over all predicted genes of each plasmid
    by_plasmid = {}
    for gid, (plasmid, n, start, stop, strand) in genes.items():
        by_plasmid.setdefault(plasmid, []).append((start, stop, strand, gid))
    operon_of = {}
    for plasmid, glist in by_plasmid.items():
        glist.sort()
        groups = []
        for i, (start, stop, strand, gid) in enumerate(glist):
            if groups:
                prev = glist[groups[-1][-1]]
                gap = start - prev[1] - 1
                if strand == prev[2] and args.min_gap <= gap <= args.max_gap:
                    groups[-1].append(i)
                    continue
            groups.append([i])
        for k, group in enumerate(groups, 1):
            letters = "".join(best[glist[i][3]]["ta_type"] if glist[i][3] in best else "x" for i in group)
            for i in group:
                operon_of[glist[i][3]] = (f"{plasmid}_op{k}", letters)

    n_count = 0
    with open(f"{args.prefix}.ta.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["plasmid", "gene_id", "start", "end", "strand", "ta_type", "best_hmm", "bitscore", "evalue",
                    "operon_id", "operon_structure", "bakta_locus", "bakta_gene", "bakta_product", "repabc", "count"])
        for gid in sorted(best, key=lambda g: (genes[g][0], genes[g][1])):
            r = best[gid]
            plasmid, _, start, stop, strand = genes[gid]
            match = bakta_match(bakta.get(plasmid, []), start, stop)
            loci = "; ".join(m[0] for m in match) or "-"
            gene = "; ".join(dict.fromkeys(m[1] for m in match if m[1])) or "-"
            product = "; ".join(dict.fromkeys(m[2] for m in match)) or "-"
            repabc = int(any(m[1] in REPABC for m in match))
            count = int(r["evalue"] < args.evalue and not repabc)
            n_count += count
            op_id, structure = operon_of[gid]
            w.writerow([plasmid, gid, start, stop, strand, r["ta_type"], r["hmm_id"], f"{r['bitscore']:.2f}",
                        f"{r['evalue']:.3g}", op_id, structure, loci, gene, product, repabc, count])

    print(f"{len(genes)} Prodigal genes in {len(seqs)} plasmids; {len(hmms)} HMMs (Z = {args.z:g}): {len(rows)} hits; "
          f"{len(best)} TA genes at E < {args.report_evalue:g}, {n_count} counted (E < {args.evalue:g}, not RepABC)")


if __name__ == "__main__":
    main()
