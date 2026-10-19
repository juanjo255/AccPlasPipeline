/*
 * Pangenome of whole genomes plus standalone accessory plasmids (Panaroo, as in
 * the publication: --clean-mode moderate --remove-invalid-genes).
 */
process PANAROO {
    tag "all"
    label 'process_panaroo'

    conda "${moduleDir}/environment.yml"

    input:
    path gffs, stageAs: 'gff/*'

    output:
    path "panaroo/gene_presence_absence.csv", emit: gpa
    path "panaroo/gene_data.csv"            , emit: gene_data
    path "panaroo/*"                        , emit: all
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    panaroo \\
        ${args} \\
        --clean-mode moderate \\
        --remove-invalid-genes \\
        --threads ${task.cpus} \\
        -i gff/* \\
        -o panaroo

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(panaroo --version 2>&1 | sed 's/^panaroo //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p panaroo
    touch panaroo/gene_presence_absence.csv panaroo/gene_data.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: stub
    END_VERSIONS
    """
}
