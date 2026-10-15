# Accessory Plasmids Identifier

Nextflow pipeline to detect accessory plasmids in *Sinorhizobium meliloti* genome assemblies.

This is result of our work done annotating and describing novel *Sinorhizobium meliloti* accessory plasmids.

## 0. Genome assembly (upstream, not part of this pipeline)

The pipeline starts from genome assemblies. We recommend assembling with
[nf-core/bacass 2.6.0](https://nf-co.re/bacass/2.6.0), for example:

```bash
nextflow run nf-core/bacass -r 2.6.0 -profile singularity \
    --input bacass_samplesheet.csv --assembly_type hybrid --outdir results_bacass
```

Assemblies from any other assembler also work.

## Input

A CSV samplesheet with one assembly per sample (see `assets/samplesheet_example.csv`):

| column     | description                                              |
|------------|----------------------------------------------------------|
| `sample`   | unique sample name, without whitespace                   |
| `assembly` | path to the assembly FASTA (`.fasta`, `.fa`, `.fna`, `.gz` accepted) |

## Usage

```bash
nextflow run accessory_plasmids -profile singularity \
    --input samplesheet.csv --genomad_db /path/to/genomad_db \
    --plsdb_fasta plsdb.fasta --plsdb_nuccore nuccore.csv --plsdb_typing typing.csv \
    --outdir results
```

Profiles: `singularity`, `apptainer`, `docker` and `conda` choose how software is provided, and `slurm` submits jobs to SLURM.
Combine them with a comma, e.g. `-profile singularity,slurm`.

## Steps

### 1. Contig length filter (`seqkit seq`)

Keeps contigs **up to `--max_contig_len` bp each** (default `1000000`). This removes the chromosome, pSymA and pSymB, leaving candidate accessory replicons.

By default, each contig ID gets the sample name appended (`contig_4` → `contig_4_TG107`),
so IDs stay unique when contigs from many genomes are pooled. Turn this off with
`--add_sample_suffix false`.

Samples with no contig under the threshold are reported with a warning and skipped in later steps.

Output: `results/01_contig_length_filter/`
- `<sample>.max_len_filtered.fasta`: retained contigs
- `<sample>.contig_lengths.tsv`: retained contig IDs and lengths

### 2. Plasmid classification (`geNomad end-to-end`)

Classifies each retained contig as plasmid, virus or chromosome with
[geNomad](https://github.com/apcamargo/genomad) 1.12.0, using default thresholds
(no `--relaxed`/`--conservative`). Score calibration is off (geNomad default), so
per-sample results are identical to running all contigs pooled.

The geNomad database is required (`--genomad_db`). The publication used **database v1.9**.
Download it once, outside the pipeline:

```bash
genomad download-database /path/to/
```

`download-database` always fetches the latest release. For reproducibility, check
`genomad_db/version.txt`; it is also recorded in `versions.yml`.

At this step the geNomad results are reported only; no contig is removed.

Output: `results/02_genomad/<sample>/`
- `<sample>_summary/<sample>_plasmid_summary.tsv`: contigs classified as plasmid (score, topology, hallmarks, conjugation and AMR genes)
- `<sample>_summary/<sample>_plasmid.fna`, `_plasmid_genes.tsv`, `_plasmid_proteins.faa`: plasmid sequences, genes and proteins
- `<sample>_summary/<sample>_virus_summary.tsv`, `_virus.fna`: contigs classified as virus
- `<sample>_aggregated_classification/<sample>_aggregated_classification.tsv`: chromosome/plasmid/virus scores for every contig
- `*.log`: geNomad logs

### 3. Related PLSDB plasmids (`select_plsdb.py`, optional `mash dist`, `seqkit grep`)

Collects public plasmids related to *S. meliloti* from [PLSDB](https://doi.org/10.1093/nar/gkae1095)
(publication: version 2024_05_31_v2). They are pooled with the sample contigs for de-replication in the next step.

A PLSDB record is selected if:
1. **Metadata (always):** `--plsdb_keyword` (default `meliloti`, case-insensitive) appears in the
   MOB-typer nearest-neighbour species (`mash_neighbor_identification`) or in the record description
   (`NUCCORE_Description`). With PLSDB 2024_05_31_v2 this gives the 224 records used in the publication.
2. **Mash (optional, `--plsdb_mash`):** a PLSDB plasmid is within `--mash_max_dist` and `--mash_max_pvalue`
   of any filtered contig (`mash dist -i` against the PLSDB sketch).

No length filter is applied to PLSDB records: megaplasmid-sized records (e.g. pSymA/pSymB) must stay in so
that clusters containing them can be discarded during de-replication.

All PLSDB inputs are files distributed with PLSDB: the sequences FASTA, `nuccore.csv` (record metadata)
and `typing.csv` (MOB-typer results). The two tables are joined on `NUCCORE_ACC`; their rows are not in the
same order. The run stops if an accession is duplicated or missing from `typing.csv`.

Output: `results/03_plsdb_related/`
- `plsdb_related.fasta`: selected PLSDB plasmids, headers trimmed to the accession
- `plsdb_related.circular.tsv`: `accession<tab>1/0` (circular/linear), used by MobMess
- `plsdb_related.sources.tsv`: why each record was selected (`metadata`, `mash` or both)
- `plsdb_related.ids.txt`: selected accessions
- `<sample>.mash_dist_plsdb.tsv`: Mash hits (only with `--plsdb_mash`)

## Parameters

| parameter             | default   | description                                    |
|-----------------------|-----------|------------------------------------------------|
| `--input`             | required  | samplesheet CSV                                |
| `--outdir`            | `results` | output directory                               |
| `--max_contig_len`    | `1000000` | maximum contig length (bp) to keep             |
| `--add_sample_suffix` | `true`    | append `_<sample>` to contig IDs               |
| `--genomad_db`        | required  | geNomad database directory                     |
| `--plsdb_fasta`       | required  | PLSDB sequences FASTA                          |
| `--plsdb_nuccore`     | required  | PLSDB `nuccore.csv`                            |
| `--plsdb_typing`      | required  | PLSDB `typing.csv` (MOB-typer)                 |
| `--plsdb_keyword`     | `meliloti`| keyword for metadata-based selection           |
| `--plsdb_mash`        | `false`   | also add PLSDB Mash hits to the filtered contigs |
| `--plsdb_sketch`      | -         | PLSDB Mash sketch (required with `--plsdb_mash`) |
| `--mash_max_dist`     | `0.1`     | `mash dist -d`                                 |
| `--mash_max_pvalue`   | `0.1`     | `mash dist -v`                                 |

## Software

- Nextflow ≥ 24.04
- seqkit 2.13.0
- geNomad 1.12.0 (database v1.9)
- Mash 2.3, Python 3.12
- PLSDB 2024_05_31_v2
