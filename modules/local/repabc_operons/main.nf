/*
 * RepABC operons per plasmid and one representative RepA (complete operon, or a
 * single RepA) for replicon typing (bin/repabc_operons.py).
 */
process REPABC_OPERONS {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path bakta_tsv
    path clusters
    path ids

    output:
    path "repabc.operons.tsv", emit: operons
    path "repabc.summary.tsv", emit: summary
    path "repabc.repA.tsv"   , emit: repa
    path "versions.yml"      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    repabc_operons.py \\
        ${args} \\
        --bakta ${bakta_tsv} \\
        --clusters ${clusters} \\
        --ids ${ids} \\
        --prefix repabc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch repabc.operons.tsv repabc.summary.tsv repabc.repA.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
