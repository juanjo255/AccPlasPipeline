/*
 * Gene families shared between accessory plasmids and chromosome / pSymA / pSymB,
 * with their COG categories (bin/gene_sharing.py).
 */
process GENE_SHARING {
    tag "all"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path gpa
    path gene_data
    path tags, stageAs: 'tags/*'
    path plasmid_ids
    path eggnog
    path gff_samples

    output:
    path "family_group_membership.csv"  , emit: membership
    path "gene_sharing_summary.tsv"     , emit: summary
    path "*.{tsv,png,svg}"              , emit: all
    path "versions.yml"                 , emit: versions

    script:
    def eggnog_arg = eggnog ? "--eggnog ${eggnog}" : ''
    """
    gene_sharing.py \\
        --gpa ${gpa} \\
        --gene-data ${gene_data} \\
        --tags tags/* \\
        --plasmids ${plasmid_ids} \\
        --gff-samples ${gff_samples} \\
        ${eggnog_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch family_group_membership.csv gene_sharing_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
