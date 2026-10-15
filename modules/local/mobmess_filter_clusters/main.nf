/*
 * Keep MobMess clusters that look like real plasmids (bin/filter_mobmess_clusters.py):
 * no fragments, all members < max_len, and (compound / maximal) >= 1 circular member;
 * maximal clusters also need >= 1 sample contig.
 */
process MOBMESS_FILTER_CLUSTERS {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path clusters
    path fasta
    path plsdb_ids
    val max_len

    output:
    path "mobmess.filtered_clusters.tsv", emit: clusters
    path "mobmess.filtered_ids.txt"     , emit: ids
    path "mobmess.filtered.fasta"       , emit: fasta
    path "versions.yml"                 , emit: versions

    script:
    """
    filter_mobmess_clusters.py \\
        --clusters ${clusters} \\
        --fasta ${fasta} \\
        --plsdb-ids ${plsdb_ids} \\
        --max-len ${max_len} \\
        --prefix mobmess

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch mobmess.filtered_clusters.tsv mobmess.filtered_ids.txt mobmess.filtered.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
