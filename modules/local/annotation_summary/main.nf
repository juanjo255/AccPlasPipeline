/*
 * Per-plasmid summary and figures of TA, defense, BGC and CAZyme annotations
 * (bin/annotation_summary.py).
 */
process ANNOTATION_SUMMARY {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path plasmid_summary
    path ta
    path defense, stageAs: 'defense/*'
    path antismash
    path cazymes
    path signalp

    output:
    path "annotation_summary.tsv", emit: summary
    path "*.tsv"                 , emit: tables
    path "*.{png,svg}"           , emit: figures, optional: true
    path "versions.yml"          , emit: versions

    script:
    def defense_files = [defense].flatten()
    def defense_arg   = defense_files ? "--defense ${defense_files.join(' ')}" : ''
    def signalp_arg   = signalp ? "--signalp ${signalp}" : ''
    """
    annotation_summary.py \\
        --plasmids ${plasmid_summary} \\
        --ta ${ta} \\
        ${defense_arg} \\
        --antismash ${antismash} \\
        --cazymes ${cazymes} \\
        ${signalp_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch annotation_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
