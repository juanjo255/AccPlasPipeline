/*
 * Toxin-antitoxin genes like the TASer web server (bin/tasmania_search.py):
 * Prodigal genes trained on all plasmids, TASmania HMMs with pyhmmer (Z = 1),
 * best HMM per gene, T/A labels and pseudo-operon structure; linked to Bakta genes.
 */
process TASMANIA {
    tag "all"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path fasta
    path bakta_tsv
    path hmm_dir
    val evalue
    val report_evalue

    output:
    path "tasmania.genes.faa", emit: genes
    path "tasmania.hits.csv", emit: hits
    path "tasmania.ta.tsv"  , emit: ta
    path "versions.yml"     , emit: versions

    script:
    """
    tasmania_search.py \\
        --fasta ${fasta} \\
        --hmm-dir ${hmm_dir} \\
        --bakta ${bakta_tsv} \\
        --evalue ${evalue} \\
        --report-evalue ${report_evalue} \\
        --threads ${task.cpus} \\
        --prefix tasmania

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: \$(python3 -c 'import pyhmmer; print(pyhmmer.__version__)')
        pyrodigal: \$(python3 -c 'import pyrodigal; print(pyrodigal.__version__)')
        tasmania_hmms: \$(ls ${hmm_dir}/*.hmm3 | wc -l)
    END_VERSIONS
    """

    stub:
    """
    touch tasmania.genes.faa tasmania.hits.csv tasmania.ta.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: stub
    END_VERSIONS
    """
}
