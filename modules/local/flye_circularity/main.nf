/*
 * Circular (1) / linear (0) flag and read depth relative to the chromosome of
 * each length-filtered contig, from Flye's assembly_info.txt (`circ.`, `cov.`).
 * IDs are renamed as in the length filter.
 */
process FLYE_CIRCULARITY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(assembly_info)
    val add_sample_suffix

    output:
    tuple val(meta), path("${meta.id}.circular.tsv"), emit: circular
    tuple val(meta), path("${meta.id}.coverage.tsv"), emit: coverage
    path "versions.yml"                             , emit: versions

    script:
    """
    flye_circularity.py \\
        --info ${assembly_info} \\
        --fasta ${fasta} \\
        --sample ${meta.id} \\
        --add-suffix ${add_sample_suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.circular.tsv ${meta.id}.coverage.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
