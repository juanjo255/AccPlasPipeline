/*
 * CAZymes with run_dbcan (CAZy DIAMOND, dbCAN HMM, dbCAN-sub), then keep those
 * supported by >= min_tools tools (bin/cazyme_filter.py). SignalP is run
 * separately (optional SIGNALP process).
 */
process DBCAN {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path proteins
    path db
    val methods
    val min_tools

    output:
    path "dbcan/overview.tsv", emit: overview
    path "dbcan/*"           , emit: all
    path "cazymes.tsv"       , emit: cazymes
    path "cazymes.faa"       , emit: faa
    path "versions.yml"      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    run_dbcan CAZyme_annotation \\
        ${args} \\
        --mode protein \\
        --input_raw_data ${proteins} \\
        --db_dir ${db} \\
        --methods ${methods} \\
        --threads ${task.cpus} \\
        --output_dir dbcan

    # The publication's DIAMOND search silently returned nothing: fail instead
    if [[ ",${methods}," == *",diamond,"* ]] && [ ! -s dbcan/diamond.out ]; then
        echo "ERROR: CAZy DIAMOND search produced no output (dbcan/diamond.out)" >&2
        exit 1
    fi

    cazyme_filter.py --overview dbcan/overview.tsv --proteins ${proteins} --min-tools ${min_tools}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(python3 -c 'import importlib.metadata as m; print(m.version("dbcan"))')
        diamond: \$(diamond version | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dbcan
    touch dbcan/overview.tsv cazymes.tsv cazymes.faa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: stub
    END_VERSIONS
    """
}
