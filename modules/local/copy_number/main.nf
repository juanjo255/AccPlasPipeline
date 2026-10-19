/*
 * Plasmid copy number: median read depth of each accessory plasmid / median depth
 * of the chromosome (samtools depth -a on reads mapped to the assembly;
 * bin/copy_number.py). pSymA/pSymB are reported as a sanity check (~1).
 */
process COPY_NUMBER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(tags), path(bam)
    path plasmid_ids
    val add_sample_suffix

    output:
    tuple val(meta), path("${meta.id}.copy_number.tsv"), emit: copy_number
    path "versions.yml"                                , emit: versions

    script:
    """
    # Stream the per-base depth (multi-GB tables are never written)
    samtools depth -a ${bam} \\
        | copy_number.py \\
            --depth - \\
            --tags ${tags} \\
            --plasmids ${plasmid_ids} \\
            --sample ${meta.id} \\
            --add-suffix ${add_sample_suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.copy_number.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}
