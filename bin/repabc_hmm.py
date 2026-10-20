#!/usr/bin/env python3
"""
RepA / RepB / RepC genes on contigs with HMMs, without a full annotation.

Genes are predicted with pyrodigal (metagenomic mode: each contig on its own) and the
proteins searched with pyhmmer against the RepABC models in assets/repabc_hmms.hmm,
using each model's gathering (GA) threshold:
  TIGR03453  partition_RepA: plasmid partitioning protein RepA   (TIGRFAMs 15.0) -> repA
  TIGR03454  partition_RepB: plasmid partitioning protein RepB   (TIGRFAMs 15.0) -> repB
  PF03428    RP-C, replication protein C N-terminal domain       (Pfam 35.0)     -> repC
  PF11800    RP-C_C, replication protein C C-terminal region     (Pfam 35.0)     -> repC
A protein hit by several models gets the gene of its best-scoring model. Proteins shorter
than --min-aa are ignored: RepABC proteins are >= ~280 aa, while ~60 aa ORFs give spurious
hits to the RepC C-terminal model.

The hits are written in Bakta's TSV layout with Bakta's product names, so that
repabc_operons.py (publication operon rules) can be run on them unchanged.

With --circular, also writes <prefix>.contig_info.tsv (contig, circular, length) for
repabc_operons.py; with --assembly-info (Flye assembly_info.txt), <prefix>.assembly_stats.tsv
with the number of contigs and of circular contigs in the whole assembly.

Outputs: <prefix>.repabc.tsv (Bakta layout), <prefix>.hits.tsv (all model hits)
"""
import argparse
import sys

import pyhmmer
import pyrodigal

GENE = {"TIGR03453": "repA", "TIGR03454": "repB", "PF03428": "repC", "PF11800": "repC"}
PRODUCT = {"repA": "plasmid partitioning protein RepA", "repB": "plasmid partitioning protein RepB",
           "repC": "plasmid replication protein RepC"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="contigs (multi-FASTA)")
    p.add_argument("--hmm", required=True, help="RepABC HMMs (assets/repabc_hmms.hmm)")
    p.add_argument("--min-aa", type=int, default=200, help="minimum protein length (aa)")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--circular", help="contig <tab> 1/0 (FLYE_CIRCULARITY)")
    p.add_argument("--assembly-info", help="Flye assembly_info.txt of the whole assembly")
    p.add_argument("--sample", help="sample name for the assembly stats table")
    p.add_argument("--prefix", default="repabc")
    return p.parse_args()


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split(None, 1)[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def main():
    args = parse_args()
    finder = pyrodigal.GeneFinder(meta=True)
    alphabet = pyhmmer.easel.Alphabet.amino()
    proteins, genes, lengths = [], {}, {}
    for contig, seq in read_fasta(args.fasta):
        lengths[contig] = len(seq)
        for i, g in enumerate(finder.find_genes(seq.encode()), 1):
            locus = f"{contig}_{i:05d}"
            genes[locus] = (contig, g.begin, g.end, "+" if g.strand == 1 else "-")
            prot = g.translate().rstrip("*")
            if len(prot) < args.min_aa:
                continue
            proteins.append(pyhmmer.easel.TextSequence(name=locus.encode(), sequence=prot).digitize(alphabet))

    with pyhmmer.plan7.HMMFile(args.hmm) as fh:
        hmms = list(fh)
    best = {}  # locus -> (score, model, evalue)
    with open(f"{args.prefix}.hits.tsv", "w") as out:
        out.write("locus\tcontig\tstart\tend\tstrand\tmodel\tscore\tevalue\n")
        for hits in pyhmmer.hmmsearch(hmms, pyhmmer.easel.DigitalSequenceBlock(alphabet, proteins),
                                      cpus=args.threads, bit_cutoffs="gathering"):
            acc = (hits.query.accession or hits.query.name)
            acc = (acc.decode() if isinstance(acc, bytes) else acc).split(".")[0]
            for h in hits.included:
                locus = h.name.decode() if isinstance(h.name, bytes) else h.name
                c, s, e, st = genes[locus]
                out.write(f"{locus}\t{c}\t{s}\t{e}\t{st}\t{acc}\t{h.score:.1f}\t{h.evalue:.2e}\n")
                if locus not in best or h.score > best[locus][0]:
                    best[locus] = (h.score, acc, h.evalue)

    rows = sorted(((genes[l][0], genes[l][1], l) for l in best))
    with open(f"{args.prefix}.repabc.tsv", "w") as out:
        out.write("#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n")
        for contig, _, locus in rows:
            c, s, e, st = genes[locus]
            score, acc, ev = best[locus]
            gene = GENE[acc]
            out.write(f"{c}\tcds\t{s}\t{e}\t{st}\t{locus}\t{gene}\t{PRODUCT[gene]}\t{acc}:{score:.1f}\n")
    if args.circular:
        with open(args.circular) as fh:
            circ = dict(l.rstrip("\n").split("\t")[:2] for l in fh if l.strip())
        with open(f"{args.prefix}.contig_info.tsv", "w") as out:
            out.write("contig\tcircular\tlength\n")
            for c, n in lengths.items():
                out.write(f"{c}\t{circ.get(c, '0')}\t{n}\n")
    if args.assembly_info:
        n_contigs = n_circ = 0
        with open(args.assembly_info) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                n_contigs += 1
                n_circ += f[3] == "Y"   # Flye column 'circ.'
        with open(f"{args.prefix}.assembly_stats.tsv", "w") as out:
            out.write("sample\tassembly_contigs\tassembly_circular_contigs\n")
            out.write(f"{args.sample or args.prefix}\t{n_contigs}\t{n_circ}\n")
    n = {g: sum(1 for l in best if GENE[best[l][1]] == g) for g in ("repA", "repB", "repC")}
    print(f"{len(proteins)} proteins on {len({v[0] for v in genes.values()})} contigs; "
          f"repA {n['repA']}, repB {n['repB']}, repC {n['repC']} genes on {len({genes[l][0] for l in best})} contigs",
          file=sys.stderr)


if __name__ == "__main__":
    main()
