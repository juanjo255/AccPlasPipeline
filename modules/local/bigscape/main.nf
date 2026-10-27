/*
 * BiG-SCAPE 2 networking of the antiSMASH regions with MIBiG reference BGCs
 * (flags of Pulliam et al. 2026: mix bin, class bins, local alignment, hybrids split, singletons kept).
 * MIBiG is passed as a reference folder of antiSMASH-processed GenBank files (BGC*.gbk), so no download
 * happens at run time.
 */
process BIGSCAPE {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path regions, stageAs: 'regions/*'
    path pfam_dir
    path mibig_dir
    val  cutoffs

    output:
    path "bigscape"              , emit: outdir
    path "bigscape/bigscape.db"  , emit: db
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Pre-create the output folder: BiG-SCAPE's own check-then-create can race on NFS
    mkdir -p bigscape

    bigscape cluster \\
        ${args} \\
        -i regions \\
        -r ${mibig_dir} \\
        --include-gbk region,BGC \\
        -p ${pfam_dir}/Pfam-A.hmm \\
        -o bigscape \\
        -db bigscape/bigscape.db \\
        -l bigscape \\
        --mix \\
        --classify class \\
        --alignment-mode local \\
        --hybrids-off \\
        --include-singletons \\
        --gcf-cutoffs ${cutoffs} \\
        -c ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: \$(python3 -c "import importlib.metadata as m; print(m.version('big-scape'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bigscape
    touch bigscape/bigscape.db
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: stub
    END_VERSIONS
    """
}
