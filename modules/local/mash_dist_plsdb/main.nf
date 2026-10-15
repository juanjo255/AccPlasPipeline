/*
 * Mash distances from each filtered contig to every PLSDB plasmid, keeping only
 * pairs within max distance (-d) and max p-value (-v). Column 1 = PLSDB accession.
 */
process MASH_DIST_PLSDB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path sketch
    val max_dist
    val max_pvalue

    output:
    tuple val(meta), path("${prefix}.mash_dist_plsdb.tsv"), emit: dist
    path "versions.yml"                                   , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # -i: compare each contig individually (the PLSDB sketch is per sequence too)
    mash dist \\
        ${args} \\
        -i \\
        -d ${max_dist} \\
        -v ${max_pvalue} \\
        -p ${task.cpus} \\
        ${sketch} \\
        ${fasta} \\
        > ${prefix}.mash_dist_plsdb.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mash_dist_plsdb.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: 2.3
    END_VERSIONS
    """
}
