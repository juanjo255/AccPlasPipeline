/*
 * Pool sample contigs and related PLSDB plasmids (sorted, IDs checked), then
 * cluster them with MobMess (MUMmer alignments, 90% identity / 90% coverage).
 */
process MOBMESS_SYSTEMS {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path fastas  , stageAs: 'fasta/*'
    path circular, stageAs: 'circular/*'

    output:
    path "mobmess_input.fasta"               , emit: fasta
    path "mobmess_input.circular.tsv"        , emit: circular
    path "${prefix}_clusters.txt"            , emit: clusters
    path "${prefix}_contigs.txt"             , emit: contigs
    path "${prefix}_systems.txt"             , emit: systems
    path "${prefix}_ani.txt"                 , emit: ani
    path "${prefix}_ani_blocks.txt"          , emit: ani_blocks
    path "${prefix}_cluster_graph.{txt,graphml}", emit: cluster_graph
    path "${prefix}_contig_graph.graphml"    , emit: contig_graph
    path "versions.yml"                      , emit: versions

    script:
    prefix   = task.ext.prefix ?: 'mobmess'
    def args = task.ext.args ?: ''
    """
    prepare_mobmess_input.py \\
        --fasta ${fastas} \\
        --circular ${circular}

    mkdir -p tmp
    mobmess systems \\
        ${args} \\
        --sequences mobmess_input.fasta \\
        --complete mobmess_input.circular.tsv \\
        --output ${prefix} \\
        --threads ${task.cpus} \\
        --tmp tmp
    # MobMess keeps --tmp (nucmer intermediates); drop it to save space in work/
    rm -rf tmp

    # MobMess and PlasX have no version number: record the installed git commits
    commit() {
        python3 -c "import json, importlib.metadata as m; print(json.loads(m.distribution('\$1').read_text('direct_url.json'))['vcs_info']['commit_id'])"
    }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobmess: \$(commit mobmess)
        plasx: \$(commit plasx)
        mummer4: \$(nucmer --version 2>&1 | tail -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'mobmess'
    """
    touch mobmess_input.fasta mobmess_input.circular.tsv
    for f in clusters.txt contigs.txt systems.txt ani.txt ani_blocks.txt \\
             cluster_graph.txt cluster_graph.graphml contig_graph.graphml; do
        touch ${prefix}_\$f
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobmess: stub
    END_VERSIONS
    """
}
