/*
 * Extract the selected PLSDB plasmids and trim headers to the accession
 * (MobMess and later steps match sequences by ID only).
 */
process PLSDB_EXTRACT {
    tag "plsdb"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fe272ab9a519cf418160471a485b5ef50ea3f571a8e4555a826f70a4d8243ae/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751' }"

    input:
    path plsdb_fasta
    path ids

    output:
    path "plsdb_related.fasta", emit: fasta
    path "versions.yml"       , emit: versions

    script:
    """
    # Exact match on sequence ID (first word of the header)
    seqkit grep --pattern-file ${ids} ${plsdb_fasta} \\
        | seqkit seq --only-id \\
        > plsdb_related.fasta

    # Every selected ID must be found in the PLSDB FASTA
    n_ids=\$(grep -c . ${ids})
    n_seq=\$(grep -c '^>' plsdb_related.fasta || true)
    if [ "\$n_ids" -ne "\$n_seq" ]; then
        echo "ERROR: \$n_ids IDs selected but \$n_seq sequences extracted from ${plsdb_fasta}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/^.*v//')
    END_VERSIONS
    """

    stub:
    """
    touch plsdb_related.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.13.0
    END_VERSIONS
    """
}
