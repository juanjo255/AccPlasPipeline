/*
 * Annotate all sequences in one Bakta run (as in the publication), keeping the
 * contig IDs so annotations can be joined back to MobMess and geNomad results.
 */
process BAKTA {
    tag "${prefix}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path fasta
    path db

    output:
    path "${prefix}.tsv"              , emit: tsv
    path "${prefix}.gff3"             , emit: gff3
    path "${prefix}.gbff"             , emit: gbff
    path "${prefix}.faa"              , emit: faa
    path "${prefix}.ffn"              , emit: ffn
    path "${prefix}.fna"              , emit: fna
    path "${prefix}.embl"             , emit: embl
    path "${prefix}.json"             , emit: json
    path "${prefix}.txt"              , emit: txt
    path "${prefix}.hypotheticals.*"  , emit: hypotheticals
    path "${prefix}.inference.tsv"    , emit: inference
    path "${prefix}.{png,svg}"        , emit: plots, optional: true
    path "versions.yml"               , emit: versions

    script:
    prefix   = task.ext.prefix ?: 'plasmids'
    def args = task.ext.args ?: ''
    """
    mkdir -p tmp
    bakta \\
        ${args} \\
        --db ${db} \\
        --keep-contig-headers \\
        --force \\
        --verbose \\
        --threads ${task.cpus} \\
        --tmp-dir tmp \\
        --prefix ${prefix} \\
        --output . \\
        ${fasta}
    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(bakta --version | sed 's/^bakta //')
        bakta_db: \$(python3 -c "import json; d=json.load(open('${db}/version.json')); print('v%s.%s %s' % (d['major'], d['minor'], d['type']))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'plasmids'
    """
    for f in tsv gff3 gbff faa ffn fna embl json txt hypotheticals.tsv hypotheticals.faa inference.tsv; do
        touch ${prefix}.\$f
    done
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: stub
    END_VERSIONS
    """
}
