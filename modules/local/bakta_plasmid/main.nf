/*
 * Bakta on one plasmid at a time (as for the publication's final annotation):
 * each plasmid's gene calls do not depend on which other plasmids are in the run.
 */
process BAKTA_PLASMID {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${meta.id}.tsv") , emit: tsv
    tuple val(meta), path("${meta.id}.faa") , emit: faa
    tuple val(meta), path("${meta.id}.gbff"), emit: gbff
    tuple val(meta), path("${meta.id}.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.{ffn,fna,embl,json,txt,hypotheticals.tsv,hypotheticals.faa,inference.tsv}"), emit: others
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p tmp
    bakta \\
        ${args} \\
        --db ${db} \\
        --keep-contig-headers \\
        --force \\
        --threads ${task.cpus} \\
        --tmp-dir tmp \\
        --prefix ${meta.id} \\
        --output . \\
        ${fasta}
    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(bakta --version | sed 's/^bakta //')
    END_VERSIONS
    """

    stub:
    """
    for f in tsv faa gbff gff3 ffn fna embl json txt hypotheticals.tsv hypotheticals.faa inference.tsv; do
        touch ${meta.id}.\$f
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: stub
    END_VERSIONS
    """
}
