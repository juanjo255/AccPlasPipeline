/*
 * BGC novelty (bin/bgc_novelty.py): known/novel from BiG-SCAPE + MIBiG, BiG-SLiCE GCFs and their
 * distance to antiSMASH-DB v4 GCFs (novelty score, antiSMASH-DB match, >= 5 BGC flag).
 */
process BGC_NOVELTY {
    tag "all"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    path bigscape_db
    path bigslice_db, stageAs: 'bigslice.db'
    path query_db   , stageAs: 'query.db'
    path asdb
    val  bigscape_cutoff
    val  match_dist
    val  min_gcf_size

    output:
    path "bgc_novelty.tsv"         , emit: bgcs
    path "gcf_novelty_bigslice.tsv", emit: gcfs
    path "versions.yml"            , emit: versions

    script:
    """
    bgc_novelty.py \\
        --bigscape-db ${bigscape_db} \\
        --bigslice-db bigslice.db \\
        --query-db query.db \\
        --asdb-db ${asdb}/result/data.db \\
        --bigscape-cutoff ${bigscape_cutoff} \\
        --match-dist ${match_dist} \\
        --min-gcf-size ${min_gcf_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch bgc_novelty.tsv gcf_novelty_bigslice.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
