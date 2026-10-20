/*
 * RepA / RepB / RepC on the filtered contigs of one sample, before MobMess: pyrodigal genes
 * + RepABC HMMs (TIGR03453, TIGR03454, PF03428, PF11800; gathering thresholds), then the
 * publication operon rules (repabc_operons.py). Also counts the contigs and circular contigs
 * of the whole assembly (Flye assembly_info.txt).
 */
process REPABC_HMM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(circular), path(assembly_info)
    path hmm

    output:
    tuple val(meta), path("${meta.id}.repabc.tsv")        , emit: genes
    tuple val(meta), path("${meta.id}.summary.tsv")       , emit: operons
    tuple val(meta), path("${meta.id}.assembly_stats.tsv"), emit: assembly_stats
    tuple val(meta), path("${meta.id}.contig_info.tsv")   , emit: contig_info
    tuple val(meta), path("${meta.id}.hits.tsv")          , emit: hits
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    repabc_hmm.py \\
        ${args} \\
        --fasta ${fasta} \\
        --hmm ${hmm} \\
        --circular ${circular} \\
        --assembly-info ${assembly_info} \\
        --sample ${meta.id} \\
        --threads ${task.cpus} \\
        --prefix ${meta.id}

    repabc_operons.py \\
        --bakta ${meta.id}.repabc.tsv \\
        --clusters ${meta.id}.contig_info.tsv \\
        --prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: \$(python3 -c 'import pyhmmer; print(pyhmmer.__version__)')
        pyrodigal: \$(python3 -c 'import pyrodigal; print(pyrodigal.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf '#Sequence Id\\tType\\tStart\\tStop\\tStrand\\tLocus Tag\\tGene\\tProduct\\tDbXrefs\\n' > ${meta.id}.repabc.tsv
    printf 'sequence\\trep_genes\\tcomplete_operon\\trepA_locus\\n' > ${meta.id}.summary.tsv
    printf 'sample\\tassembly_contigs\\tassembly_circular_contigs\\n${meta.id}\\t0\\t0\\n' > ${meta.id}.assembly_stats.tsv
    printf 'contig\\tcircular\\tlength\\n' > ${meta.id}.contig_info.tsv
    touch ${meta.id}.hits.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyhmmer: stub
    END_VERSIONS
    """
}
