/*
 * All-vs-all ANI and alignment fraction between accessory plasmids with skani
 * (one FASTA per plasmid, as in the publication), plus the AF vs ANI hexbin
 * (bin/plot_skani.py).
 */
process SKANI_ANI {
    tag "all"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path fasta

    output:
    path "skani_dist.tsv" , emit: raw
    path "plasmid_ani.tsv", emit: ani
    path "ani_vs_af.{png,svg}", emit: figures
    path "versions.yml"   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    plot_skani.py split --fasta ${fasta}

    skani dist \\
        ${args} \\
        -t ${task.cpus} \\
        --ql plasmids.list \\
        --rl plasmids.list \\
        -o skani_dist.tsv

    plot_skani.py plot --skani skani_dist.tsv
    rm -rf plasmids

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani --version | sed 's/^skani //')
    END_VERSIONS
    """

    stub:
    """
    touch skani_dist.tsv plasmid_ani.tsv ani_vs_af.png ani_vs_af.svg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: stub
    END_VERSIONS
    """
}
