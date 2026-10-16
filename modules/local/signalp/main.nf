/*
 * SignalP 6 (fast, organism "other") on the CAZyme proteins. SignalP is licensed
 * and not available from conda: it is loaded as an environment module
 * (params.signalp_module) and needs its model directory (params.signalp_model_dir).
 */
process SIGNALP {
    tag "cazymes"
    label 'process_medium'

    module params.signalp_module

    input:
    path faa
    val model_dir

    output:
    path "signalp/prediction_results.txt", emit: predictions
    path "signalp/*"                     , emit: all
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    signalp6 \\
        ${args} \\
        --fastafile ${faa} \\
        --organism other \\
        --mode fast \\
        --format txt \\
        --model_dir ${model_dir} \\
        --output_dir signalp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(head -1 signalp/prediction_results.txt | sed 's/^# SignalP-\\([0-9.]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p signalp
    touch signalp/prediction_results.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: stub
    END_VERSIONS
    """
}
