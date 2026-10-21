/*
 * Compound plasmids vs their parent plasmid backbones (MobMess Step 4 run, whose cluster IDs are
 * the ones in the publication): cointegrates (two complementary parents) and the extra sequence
 * not present in the parent, classified from RepABC HMM, MOB-typer, geNomad and Bakta hints as
 * Rep / Rep and Integrase / phage-plasmid / ICE-like / IME-like / MOB-MPF (others unclassified)
 * (bin/cointegrates.py).
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
    path repabc       , stageAs: 'repabc/*'
    path mob_biomarker
    path genomad_genes, stageAs: 'genomad_genes/*'
    path genomad_virus, stageAs: 'genomad_virus/*'

    output:
    path "cointegrates.tsv"              , emit: cointegrates
    path "compound_parents.tsv"          , emit: parents
    path "extra_segments.tsv"            , emit: segments
    path "extra_segment_genes.tsv"       , emit: segment_genes
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
        --repabc repabc/* \\
        --mob-biomarker ${mob_biomarker} \\
        --genomad-genes genomad_genes/* \\
        --genomad-virus genomad_virus/* \\
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
    touch cointegrates.tsv compound_parents.tsv extra_segments.tsv extra_segment_genes.tsv integration_summary.tsv \\
        cointegrate_junction_genes.tsv integration_map.png integration_map.svg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
