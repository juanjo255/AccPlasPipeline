/*
 * COG content of accessory plasmids: category counts, gene-content ordination
 * (Jaccard on eggNOG OGs) and COG enrichment per MobMess cluster (bin/cog_analysis.py).
 */
process COG_ANALYSIS {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path annotations
    path plasmid_summary
    val min_cluster

    output:
    path "*.tsv"       , emit: tables
    path "cog_enrichment_significant.tsv", emit: significant
    path "*.{png,svg}" , emit: figures, optional: true
    path "versions.yml", emit: versions

    script:
    """
    cog_analysis.py \\
        --annotations ${annotations} \\
        --plasmids ${plasmid_summary} \\
        --min-cluster ${min_cluster}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scikit-learn: \$(python3 -c 'import sklearn; print(sklearn.__version__)')
        statsmodels: \$(python3 -c 'import statsmodels; print(statsmodels.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch cog_counts.tsv cog_enrichment_significant.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
