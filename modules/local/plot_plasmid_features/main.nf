/*
 * Summary table and figures of the accessory plasmids, coloured by the
 * samplesheet `dataset` column (bin/plot_plasmid_features.py).
 */
process PLOT_PLASMID_FEATURES {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path features
    path rep_proteins
    path mob_type
    path repabc_summary
    path replicon_groups
    path contig_datasets

    output:
    path "plasmid_summary.tsv", emit: summary
    path "*.png"              , emit: png
    path "*.svg"              , emit: svg
    path "versions.yml"       , emit: versions

    script:
    """
    plot_plasmid_features.py \\
        --features ${features} \\
        --rep-proteins ${rep_proteins} \\
        --mobtyper ${mob_type} \\
        --repabc ${repabc_summary} \\
        --replicon-groups ${replicon_groups} \\
        --datasets ${contig_datasets}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
        seaborn: \$(python3 -c 'import seaborn; print(seaborn.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch plasmid_summary.tsv features_stackbar.png features_stackbar.svg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
