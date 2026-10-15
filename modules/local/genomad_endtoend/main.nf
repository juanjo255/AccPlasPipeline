/*
 * Classify contigs as plasmid / virus / chromosome with geNomad end-to-end
 * (default thresholds, as in the publication).
 *
 * Score calibration is left off (geNomad default), so each contig's score does
 * not depend on the other contigs in the input: running per sample gives the
 * same per-contig results as running all samples pooled.
 */
process GENOMAD_ENDTOEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}_summary/${prefix}_plasmid_summary.tsv"), emit: plasmid_summary
    tuple val(meta), path("${prefix}_summary/${prefix}_plasmid.fna")        , emit: plasmid_fasta
    tuple val(meta), path("${prefix}_summary/${prefix}_plasmid_genes.tsv")  , emit: plasmid_genes
    tuple val(meta), path("${prefix}_summary/${prefix}_plasmid_proteins.faa"), emit: plasmid_proteins
    tuple val(meta), path("${prefix}_summary/${prefix}_virus_summary.tsv")  , emit: virus_summary
    tuple val(meta), path("${prefix}_summary/${prefix}_virus.fna")          , emit: virus_fasta
    tuple val(meta), path("${prefix}_aggregated_classification/${prefix}_aggregated_classification.tsv"), emit: classification
    tuple val(meta), path("${prefix}_summary/${prefix}_summary.json")       , emit: json
    tuple val(meta), path("*.log")                                          , emit: log
    path "versions.yml"                                                     , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # geNomad names outputs after the input file stem: link it to a fixed name
    ln -s ${fasta} ${prefix}.fna

    genomad end-to-end \\
        ${args} \\
        --cleanup \\
        --threads ${task.cpus} \\
        ${prefix}.fna \\
        . \\
        ${db}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(genomad --version | sed 's/^.*version //')
        genomad_db: \$(cat ${db}/version.txt)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_summary ${prefix}_aggregated_classification
    for f in plasmid_summary.tsv plasmid.fna plasmid_genes.tsv plasmid_proteins.faa \\
             virus_summary.tsv virus.fna summary.json; do
        touch ${prefix}_summary/${prefix}_\$f
    done
    touch ${prefix}_aggregated_classification/${prefix}_aggregated_classification.tsv
    touch ${prefix}_summary.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: 1.12.0
        genomad_db: stub
    END_VERSIONS
    """
}
