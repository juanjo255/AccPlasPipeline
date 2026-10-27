/*
 * BiG-SLiCE 2 GCFs of the antiSMASH regions (threshold 0.65, as Pulliam et al. 2026).
 * Input layout built here: datasets.tsv + plasmids/<record>/<record>.regionNNN.gbk (no taxonomy).
 */
process BIGSLICE_CLUSTER {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path regions, stageAs: 'regions/*'
    path models
    val  threshold

    output:
    path "bigslice"              , emit: outdir
    path "bigslice/result/data.db", emit: db
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # One "genome" folder per plasmid record; the taxonomy column is left empty (skipped by BiG-SLiCE)
    mkdir -p input/plasmids
    printf "#dataset\\tpath\\ttaxonomy\\tdescription\\nplasmids\\tplasmids/\\t\\taccessory plasmid antiSMASH regions\\n" > input/datasets.tsv
    for f in regions/*.gbk; do
        b=\$(basename "\$f")
        rec=\${b%.region*.gbk}
        mkdir -p "input/plasmids/\$rec"
        ln -s "\$(readlink -f "\$f")" "input/plasmids/\$rec/\$b"
    done

    # Node-local temp space: on NFS, leftover .nfs* files make Python's temp-folder cleanup fail
    export TMPDIR=\$(mktemp -d "\${TMPDIR:-/tmp}/claude_bigslice.XXXXXX")
    bigslice \\
        ${args} \\
        -i input \\
        --threshold ${threshold} \\
        --program_db_folder ${models} \\
        -t ${task.cpus} \\
        bigslice
    bigslice --export-tsv bigslice/tsv --program_db_folder ${models} bigslice
    # BiG-SLiCE keeps feature-extraction intermediates; drop them to save space in work/
    rm -rf bigslice/result/tmp "\$TMPDIR"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigslice: \$(python3 -c "import importlib.metadata as m; print(m.version('bigslice'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bigslice/result
    touch bigslice/result/data.db
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigslice: stub
    END_VERSIONS
    """
}
