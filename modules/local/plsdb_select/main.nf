/*
 * Select PLSDB plasmids related to the target species (see bin/select_plsdb.py):
 * keyword in MOB-typer nearest-neighbour species or in the record description,
 * plus, optionally, all mash hits to the filtered contigs. Runs once per pipeline run.
 */
process PLSDB_SELECT {
    tag "${keyword}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path nuccore
    path typing
    val keyword
    path mash_hits, stageAs: 'mash/*'

    output:
    path "plsdb_related.ids.txt"    , emit: ids
    path "plsdb_related.sources.tsv", emit: sources
    path "plsdb_related.circular.tsv", emit: circular
    path "versions.yml"              , emit: versions

    script:
    // mash_hits is empty when --plsdb_mash is off; a single file is not a list
    def mash_files = [mash_hits].flatten()
    def mash_arg   = mash_files ? "--mash ${mash_files.join(' ')}" : ''
    """
    select_plsdb.py \\
        --nuccore ${nuccore} \\
        --typing ${typing} \\
        --keyword '${keyword}' \\
        ${mash_arg} \\
        --prefix plsdb_related

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch plsdb_related.ids.txt plsdb_related.sources.tsv plsdb_related.circular.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
