/*
 * Merge per-plasmid Bakta outputs: features table, proteins with <locus>__<plasmid>
 * headers, and GenBank records (bin/merge_plasmid_annotations.py).
 */
process MERGE_ANNOTATIONS {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path tsv , stageAs: 'tsv/*'
    path faa , stageAs: 'faa/*'
    path gbff, stageAs: 'gbff/*'

    output:
    path "plasmids.tsv"         , emit: tsv
    path "plasmids_proteins.faa", emit: faa
    path "plasmids.gbff"        , emit: gbff
    path "versions.yml"         , emit: versions

    script:
    """
    merge_plasmid_annotations.py --tsv tsv/* --faa faa/* --gbff gbff/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch plasmids.tsv plasmids_proteins.faa plasmids.gbff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
