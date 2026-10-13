# accessory_plasmids

Nextflow pipeline to detect accessory plasmids in *Sinorhizobium meliloti* genome assemblies.

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
nextflow run accessory_plasmids -profile singularity --input samplesheet.csv --outdir results
```

Profiles: `singularity`, `apptainer`, `docker` and `conda` choose how software is provided, and `slurm` submits jobs to SLURM.
Combine them with a comma, e.g. `-profile singularity,slurm`.

## Steps

### 1. Contig length filter (`seqkit seq`)

Keeps contigs **up to `--max_contig_len` bp each** (default `1000000`). This removes the
chromosome, pSymA and pSymB, leaving candidate accessory replicons. The threshold is applied per contig, not to the total genome size.

By default, each contig ID gets the sample name appended (`contig_4` → `contig_4_TG107`),
so IDs stay unique when contigs from many genomes are pooled. Turn this off with
`--add_sample_suffix false`.

Samples with no contig under the threshold are reported with a warning and skipped in later steps.

Output: `results/01_contig_length_filter/`
- `<sample>.max_len_filtered.fasta`: retained contigs
- `<sample>.contig_lengths.tsv`: retained contig IDs and lengths

## Parameters

| parameter             | default   | description                                    |
|-----------------------|-----------|------------------------------------------------|
| `--input`             | required  | samplesheet CSV                                |
| `--outdir`            | `results` | output directory                               |
| `--max_contig_len`    | `1000000` | maximum contig length (bp) to keep             |
| `--add_sample_suffix` | `true`    | append `_<sample>` to contig IDs               |

## Software

- Nextflow ≥ 24.04
- seqkit 2.13.0
