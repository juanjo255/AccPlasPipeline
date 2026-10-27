/*
 * BiG-SLiCE 2 query of the antiSMASH regions against the precomputed GCF models of the antiSMASH
 * database v4.0b2 (Zenodo 10.5281/zenodo.17969757, as_db_bigslice_out.tar.gz).
 * BiG-SLiCE writes the query report into the reference folder, so the reference result/ is linked
 * into a local folder: the shared database is read only and the report stays in work/.
 */
process BIGSLICE_QUERY {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path regions, stageAs: 'regions/*'
    path asdb
    path models

    output:
    path "query/reports"                , emit: reports
    path "query/reports/1/data.db"      , emit: db
    path "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p query
    ln -s "\$(readlink -f ${asdb})/result" query/result

    # Node-local temp space: on NFS, leftover .nfs* files make Python's temp-folder cleanup fail
    export TMPDIR=\$(mktemp -d "\${TMPDIR:-/tmp}/claude_bigslice.XXXXXX")
    bigslice \\
        ${args} \\
        --query regions \\
        --query_name plasmids \\
        --n_ranks 3 \\
        --program_db_folder ${models} \\
        -t ${task.cpus} \\
        query \\
        || echo "bigslice exited non-zero; checking the query report" >&2
    rm -rf "\$TMPDIR"
    # The report is copied before BiG-SLiCE cleans up its temp folder, which can fail after a
    # successful query: judge success by the report itself
    [ -s query/reports/1/data.db ] || { echo "ERROR: no query report (query/reports/1/data.db)" >&2; exit 1; }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigslice: \$(python3 -c "import importlib.metadata as m; print(m.version('bigslice'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p query/reports/1
    touch query/reports/1/data.db
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigslice: stub
    END_VERSIONS
    """
}
