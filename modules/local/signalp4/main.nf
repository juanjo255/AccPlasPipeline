/*
 * SignalP 4.1 (gram-, no TM model) on all accessory-plasmid proteins, as used for the
 * publication's signal-peptide analysis (Fig. 9D). SignalP is licensed and not available
 * from conda: it is loaded as an environment module (params.signalp4_module, e.g. the
 * 'spada/1.0' module that provides signalp4 at Stowers).
 */
process SIGNALP4 {
    tag "plasmids"
    label 'process_medium'

    module params.signalp4_module

    input:
    path proteins

    output:
    path "signalp4.summary"     , emit: summary
    path "signalp4.mature.fasta", emit: mature
    path "signalp4.short"       , emit: short
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    signalp4 \\
        ${args} \\
        -s notm \\
        -t gram- \\
        -f summary \\
        -n signalp4.summary \\
        -m signalp4.mature.fasta \\
        ${proteins} > signalp4.short

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(grep -m1 -o 'SignalP-[0-9.]*[a-z]*' signalp4.short | sed 's/SignalP-//')
    END_VERSIONS
    """

    stub:
    """
    touch signalp4.summary signalp4.mature.fasta signalp4.short
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: stub
    END_VERSIONS
    """
}
