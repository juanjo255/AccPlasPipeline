/*
 * Replicon groups: cluster representative RepA proteins with reference RepA of
 * Cavassim et al. 2020 (MMseqs2, >= 90% identity, >= 70% coverage) and name the
 * clusters (Rh<n> after a reference, Si<n> otherwise; bin/replicon_groups.py).
 */
process REPLICON_GROUPS {
    tag "all"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path repa
    path faa
    path references
    path names

    output:
    path "repA.faa"                 , emit: faa
    path "repA_clust_cluster.tsv"   , emit: clusters
    path "repA_clust_rep_seq.fasta" , emit: rep_seq
    path "repA.replicon_groups.tsv" , emit: groups
    path "versions.yml"             , emit: versions

    script:
    def names_arg = names ? "--names ${names}" : ''
    """
    replicon_groups.py extract \\
        --repa ${repa} \\
        --faa ${faa} \\
        --references ${references} \\
        --prefix repA

    mmseqs easy-cluster \\
        -c 0.7 \\
        --min-seq-id 0.9 \\
        --threads ${task.cpus} \\
        repA.faa \\
        repA_clust \\
        tmp
    rm -rf tmp

    replicon_groups.py name \\
        --clusters repA_clust_cluster.tsv \\
        --map repA.map.tsv \\
        --references ${references} \\
        ${names_arg} \\
        --prefix repA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: \$(mmseqs version)
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch repA.faa repA_clust_cluster.tsv repA_clust_rep_seq.fasta repA.replicon_groups.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs2: stub
    END_VERSIONS
    """
}
