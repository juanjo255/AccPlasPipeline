#!/usr/bin/env python3
"""
Keep sequences that carry plasmid markers (rules of the publication).

A sequence is kept when it is not in the exclude list and:
  1. Bakta annotates a replication or partition protein (product matches
     REP_PARTITION_REGEX, case-insensitive), or
  2. geNomad reports a MOB relaxase among its conjugation genes (rescues
     plasmids without an annotated rep protein, e.g. a MOBV plasmid made
     mostly of hypothetical proteins).

MOB-typer results are reported for every sequence but are not used to decide.
PLSDB sequences are not run through geNomad, so for them only rule 1 applies.

Review flag (does not change the decision): a kept sample sequence shorter than
--flag-max-len with read depth below --flag-min-rel-cov times its chromosome's
depth is flagged `low_coverage_small`. Small plasmids are usually high-copy, so
low depth suggests contamination, plasmid loss or a mixed culture.

Outputs:
  <prefix>.features.tsv      one row per sequence: markers, decision and reason
  <prefix>.rep_proteins.tsv  Bakta features matching the rep/partition regex
  <prefix>.ids.txt           kept sequence IDs
  <prefix>.fasta             kept sequences
  <prefix>.circular.tsv      kept sequence <tab> 1/0 (circular), for MobMess
"""
import argparse
import csv
import re
import sys

# Same pattern as the publication notebook (plasmid_feature.ipynb)
REP_PARTITION_REGEX = re.compile(
    r"(?:^(?:(?:plasmid[\s/_-]*)?replication)|partition|para\s|parb\s)", re.IGNORECASE
)
BAKTA_COLUMNS = ["Sequence Id", "Type", "Start", "Stop", "Strand", "Locus Tag", "Gene", "Product", "DbXrefs"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True, help="sequences to classify (MobMess filtered FASTA)")
    p.add_argument("--clusters", required=True, help="MobMess filtered clusters table")
    p.add_argument("--bakta", required=True, help="Bakta <prefix>.tsv")
    p.add_argument("--mobtyper", required=True, help="MOB-typer <prefix>_mob_type.tsv (--multi)")
    p.add_argument("--genomad", nargs="*", default=[], help="geNomad <sample>_plasmid_summary.tsv files")
    p.add_argument("--exclude", help="IDs to exclude after manual curation (one per line, '#' comments allowed)")
    p.add_argument("--coverage", nargs="*", default=[], help="<sample>.coverage.tsv files (Flye depth)")
    p.add_argument("--flag-max-len", type=int, default=50000, help="flag only sequences shorter than this (bp)")
    p.add_argument("--flag-min-rel-cov", type=float, default=0.5, help="flag when depth / chromosome depth is below this")
    p.add_argument("--prefix", default="plasmid_features")
    return p.parse_args()


def fasta_ids(path):
    with open(path) as fh:
        return [line[1:].split(None, 1)[0] for line in fh if line.startswith(">")]


def read_exclude(path):
    if not path:
        return set()
    ids = set()
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                ids.add(line.split()[0])
    return ids


def main():
    args = parse_args()
    seq_ids = fasta_ids(args.fasta)
    if len(set(seq_ids)) != len(seq_ids):
        sys.exit(f"ERROR: duplicated IDs in {args.fasta}")
    known = set(seq_ids)

    exclude = read_exclude(args.exclude)
    unknown = exclude - known
    if unknown:
        print(f"WARNING: excluded IDs not among the sequences: {', '.join(sorted(unknown))}", file=sys.stderr)

    with open(args.clusters, newline="") as fh:
        clusters = {r["contig"]: r for r in csv.DictReader(fh, delimiter="\t")}
    missing = known - set(clusters)
    if missing:
        sys.exit(f"ERROR: sequences missing from {args.clusters}: {', '.join(sorted(missing)[:5])}")

    # Bakta: replication / partition features (comment lines, including the header, start with '#')
    rep_rows = []
    with open(args.bakta) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            row = dict(zip(BAKTA_COLUMNS, line.rstrip("\n").split("\t")))
            if REP_PARTITION_REGEX.search(row.get("Product", "")):
                rep_rows.append(row)
    rep_by_seq = {}
    for row in rep_rows:
        label = row["Gene"] or row["Product"]
        rep_by_seq.setdefault(row["Sequence Id"], []).append(label)

    # geNomad plasmid summaries (sample contigs only)
    genomad = {}
    for path in args.genomad:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                genomad[row["seq_name"]] = row

    with open(args.mobtyper, newline="") as fh:
        mobtyper = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    coverage = {}
    for path in args.coverage:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                coverage[row["contig"]] = row

    out_cols = ["sequence", "source", "cluster", "cluster_type", "circular", "length",
                "rep_partition_proteins", "genomad_plasmid_score", "genomad_mob",
                "mobtyper_rep_type", "mobtyper_relaxase", "mobtyper_mpf", "mobtyper_orit",
                "mobtyper_mobility", "depth", "rel_coverage", "keep", "reason", "review_flag"]
    kept, flagged = [], []
    with open(f"{args.prefix}.features.tsv", "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(out_cols)
        for seq in seq_ids:
            cl = clusters[seq]
            is_plsdb = cl["source"] == "plsdb"
            reps = rep_by_seq.get(seq, [])
            g = None if is_plsdb else genomad.get(seq)  # PLSDB sequences are not run through geNomad
            mob_genes = [x for x in (g["conjugation_genes"].split(";") if g and g["conjugation_genes"] not in ("", "NA") else []) if x.startswith("MOB")]
            m = mobtyper.get(seq, {})

            if seq in exclude:
                keep, reason = 0, "excluded"
            elif reps:
                keep, reason = 1, "rep_partition"
            elif mob_genes:
                keep, reason = 1, "genomad_mob"
            else:
                keep, reason = 0, "no_marker"
            if keep:
                kept.append(seq)

            cov = coverage.get(seq)
            flag = "-"
            if keep and cov and int(cl["length"]) < args.flag_max_len and float(cov["rel_coverage"]) < args.flag_min_rel_cov:
                flag = "low_coverage_small"
                flagged.append(seq)

            writer.writerow([
                seq, cl["source"], cl["cluster"], cl["cluster_type"], cl["circular"], cl["length"],
                ";".join(reps) or "-",
                "not_tested" if is_plsdb else (g["plasmid_score"] if g else "not_plasmid"),
                ";".join(sorted(set(mob_genes))) or "-",
                m.get("rep_type(s)", "NA"), m.get("relaxase_type(s)", "NA"), m.get("mpf_type", "NA"),
                m.get("orit_type(s)", "NA"), m.get("predicted_mobility", "NA"),
                cov["depth"] if cov else "NA", cov["rel_coverage"] if cov else "NA",
                keep, reason, flag,
            ])

    with open(f"{args.prefix}.rep_proteins.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=BAKTA_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rep_rows)

    keep_set = set(kept)
    with open(f"{args.prefix}.ids.txt", "w") as out:
        out.writelines(f"{s}\n" for s in kept)
    with open(f"{args.prefix}.circular.tsv", "w") as out:
        out.writelines(f"{s}\t{clusters[s]['circular']}\n" for s in kept)
    with open(args.fasta) as fh, open(f"{args.prefix}.fasta", "w") as out:
        write = False
        for line in fh:
            if line.startswith(">"):
                write = line[1:].split(None, 1)[0] in keep_set
            if write:
                out.write(line)

    n_rescued = sum(1 for s in seq_ids if s in keep_set and s not in rep_by_seq)
    print(f"Kept {len(kept)} of {len(seq_ids)} sequences ({n_rescued} rescued by geNomad MOB; "
          f"{len(exclude & known)} excluded manually); {len(flagged)} flagged for review: {', '.join(flagged) or '-'}")


if __name__ == "__main__":
    main()
