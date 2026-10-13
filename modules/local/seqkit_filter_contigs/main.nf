/*
 * Keep only contigs shorter than or equal to a maximum length (per-contig, not
 * per-genome) and, optionally, append the sample ID to every contig ID so IDs
 * are unique when contigs from many assemblies are pooled downstream
 * (e.g. contig_4 -> contig_4_TG107).
 */
process SEQKIT_FILTER_CONTIGS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Same seqkit 2.13.0 images as the nf-core seqkit/seq module
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751' }"

    input:
    tuple val(meta), path(assembly)
    val max_len
    val add_sample_suffix

    output:
    tuple val(meta), path("${prefix}.max_len_filtered.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.contig_lengths.tsv")    , emit: lengths
    path "versions.yml"                                       , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def rename = add_sample_suffix
        ? "| seqkit replace --pattern '\$' --replacement '_${meta.id}'"
        : ''
    """
    seqkit seq \\
        --max-len ${max_len} \\
        ${assembly} \\
        ${rename} \\
        > ${prefix}.max_len_filtered.fasta

    # Table of retained contigs: id <tab> length
    seqkit fx2tab --name --only-id --length ${prefix}.max_len_filtered.fasta \\
        > ${prefix}.contig_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^.*v//')
    END_VERSIONS
    """
}
