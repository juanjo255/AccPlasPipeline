/*
 * Accessory plasmids that contain another plasmid backbone (MobMess containment, Step 4
 * run, whose cluster IDs are the ones in the publication): cointegrates (two complementary
 * plasmid backbones) and the extra segments outside the plasmid partner(s), characterized
 * by gene content, mobile genes at their ends and best match (bin/cointegrates.py).
 */
process COINTEGRATES {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path clusters
    path systems
    path ani_blocks
    path plasmid_ids
    path annotation

    output:
    path "cointegrates.tsv"              , emit: cointegrates
    path "plasmid_partners.tsv"          , emit: partners
    path "extra_segments.tsv"            , emit: segments
    path "integration_summary.tsv"       , emit: summary
    path "cointegrate_junction_genes.tsv", emit: genes
    path "integration_map.{png,svg}"     , emit: figures
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    cointegrates.py \\
        --clusters ${clusters} \\
        --systems ${systems} \\
        --blocks ${ani_blocks} \\
        --plasmids ${plasmid_ids} \\
        --annotation ${annotation} \\
        --min-partner-len ${params.partner_min_len} \\
        --min-part ${params.cointegrate_min_part} \\
        --max-overlap ${params.cointegrate_max_overlap} \\
        --min-union ${params.cointegrate_min_union} \\
        --min-segment ${params.segment_min_bp} \\
        --window ${params.junction_window} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch cointegrates.tsv plasmid_partners.tsv extra_segments.tsv integration_summary.tsv \\
        cointegrate_junction_genes.tsv integration_map.png integration_map.svg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
