#!/usr/bin/env python3
"""
Clean genome of one strain: main replicons + retained accessory plasmids (+ unplaced pieces),
with NCBI source modifiers in the FASTA headers.

`select` decides the role of every contig of the assembly:
  chromosome        replicon tag Chromosome                  [location=chromosome]
  pSymA / pSymB     replicon tags pSymA / pSymB              [plasmid-name=pSymA] / [plasmid-name=pSymB]
  accessory plasmid in the final accessory plasmid set      [plasmid-name=p<strain><a,b,...>], letters by
                                                            decreasing length
  unplaced          <replicon>_fragment tags and contigs of the Step 4b review table (genuine pieces,
                    kept without location / plasmid name)
  removed           everything else (debris, IS circles, contigs not retained), written to
                    <sample>.removed.fasta with the reason
An optional curation table (sample, contig, action keep|drop, name, note) overrides any decision;
`keep` with a name starting with 'p' makes it a named plasmid, `keep` without a name an unplaced contig.
Topology comes from Flye (circ. = Y); [completeness=complete] is given to circular chromosome/plasmids.

Writes <sample>.clean.plain.fasta (IDs only, for dnaapler), <sample>.linear.txt (contigs dnaapler
must not rotate), <sample>.clean.tsv (one row per contig) and <sample>.removed.fasta.

`rotate` re-orients the circular replicons: the chromosome comes from dnaapler (dnaA, given as
--chromosome); circular plasmids (pSymA, pSymB, accessory) start at the RepA of their RepABC operon
found by the RepABC HMMs (repabc_hmm.py + repabc_operons.py on the clean sequences; the RepA of a
complete operon, else the best-scoring RepA), reverse-complemented when RepA is on the minus strand.
dnaapler's own RepA database does not cover the ParA-like RepA of RepABC replicons. Linear contigs,
unplaced contigs and plasmids without RepA keep their assembled orientation. Writes
<sample>.reoriented.fasta and <sample>.rotation.tsv.

`finalize` adds the NCBI source modifiers to the (re-oriented) sequences: <sample>.clean.fasta, and
writes <sample>.summary.tsv (one row: replicons found, topology, sizes, notes).
"""
import argparse
import csv
import sys


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


def write_fasta(fh, header, seq):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def letters(i):
    """0 -> a, 25 -> z, 26 -> aa ..."""
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(97 + r) + s
    return s


def select(args):
    sfx = f"_{args.sample}"
    rename = (lambda c: c if (not args.add_suffix or c.endswith(sfx)) else c + sfx)
    seqs = {rename(n): s for n, s in read_fasta(args.assembly)}
    circular = {}
    with open(args.assembly_info) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            circular[rename(f[0])] = f[3] == "Y"
    tags = {}
    with open(args.tags) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tags[r["contig"]] = r["tag"]
    with open(args.plasmids) as fh:
        final = {l.strip() for l in fh if l.strip()}
    review = {}
    if args.review:
        with open(args.review) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                review[r["contig"]] = r.get("review_reason", "review")
    curation = {}
    if args.curation:
        with open(args.curation) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["sample"] == args.sample:
                    curation[r["contig"]] = r

    rows = []
    for c, seq in seqs.items():
        tag = tags.get(c, "Other")
        role, name, reason = "removed", "", "not a main replicon nor a retained accessory plasmid"
        if tag == "Chromosome":
            role, name, reason = "chromosome", "chromosome", "replicon tag"
        elif tag in ("pSymA", "pSymB"):
            role, name, reason = "plasmid", tag, "replicon tag"
        elif c in final:
            role, reason = "accessory", "final accessory plasmid set"
        elif tag.endswith("_fragment"):
            role, reason = "unplaced", f"replicon piece ({tag})"
        elif c in review:
            role, reason = "unplaced", f"Step 4b review: {review[c]}"
        cur = curation.get(c)
        if cur:
            act = cur.get("action", "").strip().lower()
            if act == "drop":
                role, name = "removed", ""
            elif act == "keep":
                nm = (cur.get("name") or "").strip()
                if nm == "chromosome":
                    role, name = "chromosome", "chromosome"
                elif nm.startswith("p"):
                    role, name = "plasmid", nm
                elif role == "removed":
                    role = "unplaced"
            reason = f"curation: {cur.get('note') or act}"
        rows.append({"contig": c, "length": len(seq), "replicon_tag": tag, "role": role, "name": name,
                     "topology": "circular" if circular.get(c) else "linear", "reason": reason})

    # Accessory plasmids: p<strain><letter>, by decreasing length
    acc = sorted((r for r in rows if r["role"] == "accessory"), key=lambda r: -r["length"])
    for i, r in enumerate(acc):
        r["name"] = f"p{args.sample}{letters(i)}"
        r["role"] = "plasmid"

    order = {"chromosome": 0, "plasmid": 1, "unplaced": 2, "removed": 3}
    rows.sort(key=lambda r: (order[r["role"]], r["name"] not in ("pSymA", "pSymB"), r["name"], -r["length"]))
    with open(f"{args.sample}.clean.plain.fasta", "w") as kept, open(f"{args.sample}.removed.fasta", "w") as rem, \
         open(f"{args.sample}.linear.txt", "w") as lin:
        for r in rows:
            if r["role"] == "removed":
                write_fasta(rem, f"{r['contig']} {r['reason']}", seqs[r["contig"]])
            else:
                write_fasta(kept, r["contig"], seqs[r["contig"]])
                if r["topology"] == "linear":
                    lin.write(r["contig"] + "\n")
    with open(f"{args.sample}.clean.tsv", "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["contig", "length", "replicon_tag", "role", "name", "topology", "reason"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"[{args.sample}] " + ", ".join(f"{k}: {sum(r['role'] == k for r in rows)}" for k in order))


def finalize(args):
    with open(args.table) as fh:
        rows = {r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")}
    seqs = dict(read_fasta(args.fasta))
    kept = [r for r in rows.values() if r["role"] != "removed"]
    missing = [r["contig"] for r in kept if r["contig"] not in seqs]
    if missing:
        sys.exit(f"ERROR: [{args.sample}] sequences missing after re-orientation: {missing}")
    with open(f"{args.sample}.clean.fasta", "w") as out:
        for r in kept:
            mods = []
            if r["role"] == "chromosome":
                mods.append("[location=chromosome]")
            elif r["role"] == "plasmid":
                mods.append(f"[plasmid-name={r['name']}]")
            mods.append(f"[topology={r['topology']}]")
            if r["role"] in ("chromosome", "plasmid") and r["topology"] == "circular":
                mods.append("[completeness=complete]")
            write_fasta(out, f"{r['contig']} " + " ".join(mods), seqs[r["contig"]])

    # Summary and notes for the submission
    get = lambda role, name=None: [r for r in kept if r["role"] == role and (name is None or r["name"] == name)]
    chrom, psa, psb = get("chromosome"), get("plasmid", "pSymA"), get("plasmid", "pSymB")
    acc = [r for r in get("plasmid") if r["name"] not in ("pSymA", "pSymB")]
    notes = []
    if not psb and chrom and max(int(r["length"]) for r in chrom) > 5_000_000:
        notes.append("no pSymB and a > 5 Mb chromosome: possible chromosome-pSymB fusion")
    for role, lst in (("chromosome", chrom), ("pSymA", psa), ("pSymB", psb)):
        if not lst:
            notes.append(f"{role} not found")
        elif any(r["topology"] == "linear" for r in lst):
            notes.append(f"{role} linear")
    if get("unplaced"):
        notes.append(f"{len(get('unplaced'))} unplaced contig(s)")
    with open(f"{args.sample}.summary.tsv", "w") as out:
        out.write("sample\tchromosome_bp\tpSymA_bp\tpSymB_bp\taccessory_plasmids\taccessory_circular\t"
                  "unplaced\tremoved\ttotal_bp\tnotes\n")
        bp = lambda lst: ",".join(r["length"] for r in lst) or "-"
        out.write("\t".join(map(str, [
            args.sample, bp(chrom), bp(psa), bp(psb), len(acc), sum(r["topology"] == "circular" for r in acc),
            len(get("unplaced")), sum(r["role"] == "removed" for r in rows.values()),
            sum(len(seqs[r["contig"]]) for r in kept), "; ".join(notes) or "-"])) + "\n")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def rotate(args):
    with open(args.table) as fh:
        rows = {r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")}
    seqs = dict(read_fasta(args.fasta))
    chrom = dict(read_fasta(args.chromosome)) if args.chromosome else {}
    genes = {}
    with open(args.repabc) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#") or l.startswith("#Sequence")), delimiter="\t"):
            if r["Gene"] == "repA":
                genes[r["Locus Tag"]] = r
    complete = set()
    if args.operons:
        with open(args.operons) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("complete_operon") == "1" and r["Gene"] == "repA":
                    complete.add(r["Locus Tag"])
    log = []
    with open(f"{args.sample}.reoriented.fasta", "w") as out:
        for c, seq in seqs.items():
            r = rows[c]
            how = "as assembled"
            if r["role"] == "chromosome" and c in chrom:
                seq, how = chrom[c], "dnaapler (dnaA)"
            elif r["role"] == "plasmid" and r["topology"] == "circular":
                cand = [g for g in genes.values() if g["#Sequence Id"] == c]
                if cand:
                    score = lambda g: float(g["DbXrefs"].split(":")[-1])
                    g = max(cand, key=lambda g: (g["Locus Tag"] in complete, score(g)))
                    start, stop, L = int(g["Start"]), int(g["Stop"]), len(seq)
                    if g["Strand"] == "+":
                        seq = seq[start - 1:] + seq[:start - 1]
                    else:
                        rc = revcomp(seq)
                        new = L - stop        # 0-based start of repA on the reverse complement
                        seq = rc[new:] + rc[:new]
                    how = (f"RepA ({'complete operon' if g['Locus Tag'] in complete else 'best hit'}) at "
                           f"{start}-{stop} ({g['Strand']})")
                else:
                    how = "no RepA: as assembled"
            write_fasta(out, c, seq)
            log.append((c, r["role"], r["name"], r["topology"], how))
    with open(f"{args.sample}.rotation.tsv", "w") as out:
        out.write("contig\trole\tname\ttopology\treorientation\n")
        for x in log:
            out.write("\t".join(x) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("select")
    a.add_argument("--assembly", required=True)
    a.add_argument("--assembly-info", required=True, help="Flye assembly_info.txt (topology)")
    a.add_argument("--tags", required=True, help="<sample>.replicon_tags.tsv")
    a.add_argument("--plasmids", required=True, help="final accessory plasmid IDs (all samples)")
    a.add_argument("--review", help="linear_marker_review.tsv (Step 4b): unplaced contigs")
    a.add_argument("--curation", help="sample, contig, action (keep|drop), name, note")
    a.add_argument("--sample", required=True)
    a.add_argument("--add-suffix", action="store_true", help="append _<sample> to contig IDs (as Step 1)")
    r = sub.add_parser("rotate")
    r.add_argument("--fasta", required=True, help="<sample>.clean.plain.fasta")
    r.add_argument("--table", required=True, help="<sample>.clean.tsv")
    r.add_argument("--repabc", required=True, help="repabc_hmm.py genes on the clean sequences")
    r.add_argument("--operons", help="repabc_operons.py <prefix>.operons.tsv")
    r.add_argument("--chromosome", help="dnaapler-reoriented chromosome FASTA")
    r.add_argument("--sample", required=True)
    b = sub.add_parser("finalize")
    b.add_argument("--fasta", required=True, help="re-oriented clean sequences (IDs only)")
    b.add_argument("--table", required=True, help="<sample>.clean.tsv")
    b.add_argument("--sample", required=True)
    args = p.parse_args()
    {"select": select, "rotate": rotate, "finalize": finalize}[args.cmd](args)


if __name__ == "__main__":
    main()
