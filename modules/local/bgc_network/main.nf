/*
 * Gephi node/edge tables of the BiG-SCAPE mix-bin network (bin/bgc_gephi_export.py), annotated with
 * BGC novelty and the plasmid of origin (final MobMess cluster of Step 8, replicon group of Step 7).
 */
process BGC_NETWORK {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path bigscape
    path bgc_novelty
    path gcf_novelty
    path plasmid_clusters
    path replicon_groups
    val  cutoff
    val  top_clusters

    output:
    path "nodes.csv"   , emit: nodes
    path "edges.csv"   , emit: edges
    path "bridges.tsv" , emit: bridges
    path "versions.yml", emit: versions

    script:
    def groups_arg = replicon_groups ? "--replicon-groups ${replicon_groups}" : ''
    """
    # BiG-SCAPE names the output folder <label>_<timestamp>_c<cutoff>
    cdir=\$(ls -d ${bigscape}/output_files/*_c${cutoff} | head -1)
    bgc_gephi_export.py \\
        --network \$cdir/mix/mix_c${cutoff}.network \\
        --clustering \$cdir/mix/mix_clustering_c${cutoff}.tsv \\
        --annotations \$cdir/record_annotations.tsv \\
        --bgc-novelty ${bgc_novelty} \\
        --gcf-novelty ${gcf_novelty} \\
        --plasmid-clusters ${plasmid_clusters} \\
        ${groups_arg} \\
        --top-clusters ${top_clusters}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch nodes.csv edges.csv bridges.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
