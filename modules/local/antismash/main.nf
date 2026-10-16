/*
 * antiSMASH on all plasmid GenBank records (genes from the Bakta annotation).
 */
process ANTISMASH {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path gbff
    path databases

    output:
    path "antismash"                         , emit: outdir
    path "antismash/*.region*.gbk"           , emit: regions, optional: true
    path "versions.yml"                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    antismash \\
        ${args} \\
        --cb-general \\
        --cc-mibig \\
        --taxon bacteria \\
        --cpus ${task.cpus} \\
        --databases ${databases} \\
        --output-dir antismash \\
        ${gbff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/^antiSMASH //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p antismash
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: stub
    END_VERSIONS
    """
}
