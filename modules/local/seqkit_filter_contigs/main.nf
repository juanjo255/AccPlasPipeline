/*
 * Keep only contigs shorter than or equal to a maximum length (per-contig, not
 * per-genome) and, optionally, append the sample ID to every contig ID so IDs
 * are unique when contigs from many assemblies are pooled downstream
 * (e.g. contig_4 -> contig_4_TG107). IDs that already end in _<sample> (e.g.
 * renamed or polished assemblies) are kept as they are. Only the first word of
 * each header is used.
 */
process SEQKIT_FILTER_CONTIGS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(assembly)
    val max_len
    val add_sample_suffix

    output:
    tuple val(meta), path("${prefix}.max_len_filtered.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.contig_lengths.tsv")    , emit: lengths
    path "versions.yml"                                       , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = add_sample_suffix ? "_${meta.id}" : ''
    """
    # Keep the first word of each header; append _<sample> unless already present
    # (plain string comparison: sample names can contain '.' or '-')
    seqkit seq \\
        --max-len ${max_len} \\
        --only-id \\
        ${assembly} \\
        | awk -v sfx='${suffix}' '
            /^>/ {
                id = substr(\$1, 2)
                if (sfx != "" && substr(id, length(id) - length(sfx) + 1) != sfx) id = id sfx
                print ">" id
                next
            }
            { print }' \\
        > ${prefix}.max_len_filtered.fasta

    # Table of retained contigs: id <tab> length
    seqkit fx2tab --name --only-id --length ${prefix}.max_len_filtered.fasta \\
        > ${prefix}.contig_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^.*v//')
    END_VERSIONS
    """
}
