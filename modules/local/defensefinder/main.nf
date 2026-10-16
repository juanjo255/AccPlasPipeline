/*
 * DefenseFinder (+ AntiDefenseFinder) on the ordered proteins of one plasmid.
 */
process DEFENSEFINDER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(faa)
    path models_dir

    output:
    tuple val(meta), path("${meta.id}_defense_finder_systems.tsv"), emit: systems
    tuple val(meta), path("${meta.id}_defense_finder_genes.tsv")  , emit: genes
    tuple val(meta), path("${meta.id}_defense_finder_hmmer.tsv")  , emit: hmmer
    path "versions.yml"                                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Output files are named after the input file: use the plasmid ID
    [ "${faa}" = "${meta.id}.faa" ] || ln -s ${faa} ${meta.id}.faa
    defense-finder run \\
        ${args} \\
        --workers ${task.cpus} \\
        --db-type ordered_replicon \\
        --coverage 0.4 \\
        --antidefensefinder \\
        --models-dir ${models_dir} \\
        --out-dir . \\
        ${meta.id}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: \$(python3 -c 'import importlib.metadata as m; print(m.version("mdmparis-defense-finder"))')
        macsyfinder: \$(python3 -c 'import importlib.metadata as m; print(m.version("macsyfinder"))')
        defense-finder-models: \$(grep '^vers' ${models_dir}/defense-finder-models/metadata.yml | sed 's/^vers: *//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_defense_finder_systems.tsv ${meta.id}_defense_finder_genes.tsv ${meta.id}_defense_finder_hmmer.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: stub
    END_VERSIONS
    """
}
