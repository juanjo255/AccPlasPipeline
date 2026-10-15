/*
 * Interactive MobMess containment network (all clusters, before filtering):
 * nodes coloured by cluster type and sized by length, annotated with GC,
 * source (sample / PLSDB), geNomad call and the cluster filter result.
 */
process MOBMESS_NETWORK {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path graph
    path fasta
    path plsdb_ids
    path filtered_ids
    path genomad_summaries, stageAs: 'genomad/*'

    output:
    path "mobmess_network.html"     , emit: html
    path "mobmess_network.graphml"  , emit: graphml
    path "mobmess_network.nodes.tsv", emit: nodes
    path "versions.yml"             , emit: versions

    script:
    // genomad_summaries is a list (one plasmid_summary.tsv per sample) or a single file
    def genomad_files = [genomad_summaries].flatten()
    def genomad_arg   = genomad_files ? "--genomad ${genomad_files.join(' ')}" : ''
    """
    plot_mobmess_network.py \\
        --graph ${graph} \\
        --fasta ${fasta} \\
        --plsdb-ids ${plsdb_ids} \\
        --filtered-ids ${filtered_ids} \\
        ${genomad_arg} \\
        --prefix mobmess_network

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        networkx: \$(python3 -c 'import networkx; print(networkx.__version__)')
        ipysigma: \$(python3 -c 'import importlib.metadata as m; print(m.version("ipysigma"))')
    END_VERSIONS
    """

    stub:
    """
    touch mobmess_network.html mobmess_network.graphml mobmess_network.nodes.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ipysigma: stub
    END_VERSIONS
    """
}
