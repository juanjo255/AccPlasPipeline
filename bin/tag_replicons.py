#!/usr/bin/env python3
"""
Tag every contig of an assembly as Chromosome, pSymA, pSymB or Other.

Publication rule (accessory_pls_pipeline.sh): the assembly is aligned to a reference
genome (S. meliloti GCF_037023865) with
    minimap2 -x asm5 --end-bonus 100 --secondary=no <reference> <assembly>
each contig keeps its longest alignment block, and each reference replicon is given
to the single contig with the longest such block. This reproduces all 346
publication strains (1,829 contigs).

Safeguard (new): the tag is kept only if the contig covers at least --min-ref-cov %
of that replicon (alignments merged on the reference). Without it, when a replicon
is fused to another one or missing, the tag slides to an unrelated contig: in the
publication two plasmids were tagged pSymB while covering < 1% of pSymB. All
other publication tags cover >= 41% of their replicon. Use --min-ref-cov 0 to
reproduce the publication exactly.
Contig IDs are renamed as in the contig length filter (`_<sample>` added if missing).

Outputs:
  <sample>.replicon_tags.tsv     contig, tag, length, reference, reference_coverage
  <sample>.main_replicons.fasta  contigs tagged Chromosome / pSymA / pSymB (renamed IDs)
"""
import argparse
import collections
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="minimap2 PAF of the assembly against the reference")
    p.add_argument("--assembly", required=True)
    p.add_argument("--ref-tags", required=True, help="reference replicon <tab> tag (e.g. NZ_CP146208.1 Chromosome)")
    p.add_argument("--sample", required=True)
    p.add_argument("--add-suffix", choices=["true", "false"], required=True)
    p.add_argument("--min-ref-cov", type=float, default=30, help="minimum coverage of the reference replicon (%%)")
    return p.parse_args()


def merged_length(intervals):
    total, cur_start, cur_end = 0, None, None
    for start, end in sorted(intervals):
        if cur_end is None or start > cur_end:
            if cur_end is not None:
                total += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    return total + (cur_end - cur_start if cur_end is not None else 0)


def main():
    args = parse_args()
    suffix = f"_{args.sample}"

    lengths, seqs, name = {}, {}, None
    with open(args.assembly) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                lengths[name] = 0
                seqs[name] = []
            elif name is not None:
                lengths[name] += len(line.strip())
                seqs[name].append(line.strip())

    ref_tags = {}
    with open(args.ref_tags) as fh:
        for line in fh:
            if line.strip():
                ref, tag = line.rstrip("\n").split("\t")[:2]
                ref_tags[ref] = tag

    rows, blocks, ref_len = [], collections.defaultdict(list), {}
    with open(args.paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, t, tstart, tend, block = f[0], f[5], int(f[7]), int(f[8]), int(f[10])
            if t not in ref_tags:
                sys.exit(f"ERROR: reference sequence {t} has no tag in {args.ref_tags}")
            rows.append((block, q, t))
            blocks[(q, t)].append((tstart, tend))
            ref_len[t] = int(f[6])

    # Longest block per contig, then each replicon to the contig with the longest block
    rows.sort(key=lambda r: -r[0])  # stable, like csvtk sort -k 11:nr
    best_of_contig = {}
    for block, q, t in rows:
        best_of_contig.setdefault(q, (block, t))
    assigned, taken = {}, set()
    for q, (block, t) in sorted(best_of_contig.items(), key=lambda kv: -kv[1][0]):
        if t not in taken:
            taken.add(t)
            assigned[q] = t

    counts = collections.Counter()
    with open(f"{args.sample}.replicon_tags.tsv", "w") as out, open(f"{args.sample}.main_replicons.fasta", "w") as fa:
        out.write("contig\ttag\tlength\treference\treference_coverage\n")
        for contig, length in lengths.items():
            ref, tag, cov = "-", "Other", 0.0
            if contig in assigned:
                ref = assigned[contig]
                cov = 100 * merged_length(blocks[(contig, ref)]) / ref_len[ref]
                if cov >= args.min_ref_cov:
                    tag = ref_tags[ref]
                else:
                    print(f"WARNING: [{args.sample}] {contig} best matches {ref_tags[ref]} but covers only "
                          f"{cov:.1f}% of it: tagged Other", file=sys.stderr)
            counts[tag] += 1
            new_id = contig if args.add_suffix == "false" or contig.endswith(suffix) else contig + suffix
            out.write(f"{new_id}\t{tag}\t{length}\t{ref}\t{cov:.1f}\n")
            if tag != "Other":
                seq = "".join(seqs[contig])
                fa.write(f">{new_id}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
    print(f"[{args.sample}] " + ", ".join(f"{t}: {counts[t]}" for t in ("Chromosome", "pSymA", "pSymB", "Other")))


if __name__ == "__main__":
    main()
