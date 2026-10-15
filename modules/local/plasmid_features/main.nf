/*
 * Keep sequences with plasmid markers (bin/plasmid_features.py): a Bakta
 * replication/partition protein, or a geNomad MOB relaxase; minus a manual
 * exclude list. MOB-typer results are reported alongside, and small plasmids
 * with low read depth (Flye) are flagged for manual review.
 */
process PLASMID_FEATURES {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path fasta
    path clusters
    path bakta_tsv
    path mob_type
    path genomad_summaries, stageAs: 'genomad/*'
    path exclude
    path coverage, stageAs: 'coverage/*'
    val flag_max_len
    val flag_min_rel_cov

    output:
    path "plasmid_features.features.tsv"    , emit: features
    path "plasmid_features.rep_proteins.tsv", emit: rep_proteins
    path "plasmid_features.ids.txt"         , emit: ids
    path "plasmid_features.fasta"           , emit: fasta
    path "versions.yml"                     , emit: versions

    script:
    def genomad_files = [genomad_summaries].flatten()
    def genomad_arg   = genomad_files ? "--genomad ${genomad_files.join(' ')}" : ''
    def exclude_arg   = exclude ? "--exclude ${exclude}" : ''
    def coverage_files = [coverage].flatten()
    def coverage_arg  = coverage_files ? "--coverage ${coverage_files.join(' ')}" : ''
    """
    plasmid_features.py \\
        --fasta ${fasta} \\
        --clusters ${clusters} \\
        --bakta ${bakta_tsv} \\
        --mobtyper ${mob_type} \\
        ${genomad_arg} \\
        ${exclude_arg} \\
        ${coverage_arg} \\
        --flag-max-len ${flag_max_len} \\
        --flag-min-rel-cov ${flag_min_rel_cov} \\
        --prefix plasmid_features

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch plasmid_features.features.tsv plasmid_features.rep_proteins.tsv plasmid_features.ids.txt plasmid_features.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
