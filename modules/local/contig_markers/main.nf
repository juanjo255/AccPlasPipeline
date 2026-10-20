/*
 * Plasmid markers of every filtered sample contig (RepABC HMMs, geNomad, circularity,
 * assembly circular-contig counts, MobMess cluster) and the linear contigs with RepABC genes
 * dropped by the MobMess filter, for user review (bin/contig_markers.py).
 */
process CONTIG_MARKERS {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path datasets
    path contig_info   , stageAs: 'info/*'
    path repabc        , stageAs: 'repabc/*'
    path operons       , stageAs: 'operons/*'
    path assembly_stats, stageAs: 'stats/*'
    path genomad       , stageAs: 'genomad/*'
    path coverage      , stageAs: 'coverage/*'
    path clusters
    path kept
    path tags          , stageAs: 'tags/*'
    path graph         , stageAs: 'graph/*'
    val genomad_min_score

    output:
    path "contig_markers.tsv"       , emit: markers
    path "linear_marker_review.tsv" , emit: review
    path "versions.yml"             , emit: versions

    script:
    def tag_arg   = tags ? "--tags tags/*" : ''
    def graph_arg = graph ? "--graph graph/*" : ''
    """
    contig_markers.py \\
        --datasets ${datasets} \\
        --contig-info info/* \\
        --repabc repabc/* \\
        --operons operons/* \\
        --assembly-stats stats/* \\
        --genomad genomad/* \\
        --coverage coverage/* \\
        --clusters ${clusters} \\
        --kept ${kept} \\
        --genomad-min-score ${genomad_min_score} \\
        ${tag_arg} \\
        ${graph_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch contig_markers.tsv linear_marker_review.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
