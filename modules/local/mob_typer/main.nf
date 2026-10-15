/*
 * MOB-typer on every sequence separately (--multi): replicon, relaxase, MPF and
 * oriT types and predicted mobility.
 */
process MOB_TYPER {
    tag "${prefix}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path fasta
    path db

    output:
    path "${prefix}_mob_type.tsv"           , emit: mob_type
    path "${prefix}_mobtyper_MGE_report.tsv", emit: mge_report
    path "${prefix}_mobtyper_biomarker.tsv" , emit: biomarker
    path "versions.yml"                     , emit: versions

    script:
    prefix   = task.ext.prefix ?: 'plasmids'
    def args = task.ext.args ?: ''
    """
    mob_typer \\
        ${args} \\
        --multi \\
        -n ${task.cpus} \\
        -d ${db} \\
        --infile ${fasta} \\
        --out_file ${prefix}_mob_type.tsv \\
        --mge_report_file ${prefix}_mobtyper_MGE_report.tsv \\
        --biomarker_report_file ${prefix}_mobtyper_biomarker.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mob_suite: \$(mob_typer --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'plasmids'
    """
    touch ${prefix}_mob_type.tsv ${prefix}_mobtyper_MGE_report.tsv ${prefix}_mobtyper_biomarker.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mob_suite: stub
    END_VERSIONS
    """
}
