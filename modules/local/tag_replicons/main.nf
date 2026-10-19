/*
 * Tag contigs as Chromosome / pSymA / pSymB / Other by alignment to a reference
 * genome (minimap2 asm5, publication rule plus a minimum replicon coverage;
 * bin/tag_replicons.py). Also writes the main replicon sequences.
 */
process TAG_REPLICONS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(assembly)
    path reference
    path reference_tags
    val add_sample_suffix
    val min_ref_cov

    output:
    tuple val(meta), path("${meta.id}.replicon_tags.tsv")   , emit: tags
    tuple val(meta), path("${meta.id}.main_replicons.fasta"), emit: main
    path "versions.yml"                                     , emit: versions

    script:
    """
    minimap2 -x asm5 --end-bonus 100 --secondary=no -t ${task.cpus} ${reference} ${assembly} > ${meta.id}.paf

    tag_replicons.py \\
        --paf ${meta.id}.paf \\
        --assembly ${assembly} \\
        --ref-tags ${reference_tags} \\
        --sample ${meta.id} \\
        --add-suffix ${add_sample_suffix} \\
        --min-ref-cov ${min_ref_cov}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.replicon_tags.tsv ${meta.id}.main_replicons.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: stub
    END_VERSIONS
    """
}
