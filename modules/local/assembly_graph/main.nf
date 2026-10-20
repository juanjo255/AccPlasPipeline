/*
 * Flye assembly-graph context of each contig (bin/assembly_graph.py): connected component,
 * other contigs in it, isolated or not, and repeat edges on its path.
 */
process ASSEMBLY_GRAPH {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gfa), path(assembly_info)
    val add_sample_suffix

    output:
    tuple val(meta), path("${meta.id}.graph.tsv"), emit: graph
    path "versions.yml"                          , emit: versions

    script:
    def suffix = add_sample_suffix ? '--add-sample-suffix' : ''
    """
    assembly_graph.py \\
        --gfa ${gfa} \\
        --assembly-info ${assembly_info} \\
        --sample ${meta.id} \\
        ${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'contig\\tgraph_component\\tcomponent_contigs\\tisolated_component\\tcomponent_bp\\tcomponent_edges\\tpath_edges\\trepeat_edges\\tmax_edge_depth_ratio\\n' > ${meta.id}.graph.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
