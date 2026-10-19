/*
 * Functional annotation (COG categories, orthologous groups) of the accessory
 * plasmid proteins with eggNOG-mapper (MMseqs2 search, as in the publication).
 */
process EGGNOG_MAPPER {
    tag "all"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    path proteins
    path db

    output:
    path "plasmids_eggnog.emapper.annotations"   , emit: annotations
    path "plasmids_eggnog.emapper.*"             , emit: all
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p tmp
    emapper.py \\
        ${args} \\
        --data_dir ${db} \\
        -i ${proteins} \\
        --itype proteins \\
        -m mmseqs \\
        --decorate_gff yes \\
        --cpu ${task.cpus} \\
        --temp_dir tmp \\
        -o plasmids_eggnog \\
        --output_dir .
    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | grep -o 'emapper-[0-9.]*' | head -1 | sed 's/emapper-//')
        eggnog_db: \$(emapper.py --data_dir ${db} --version 2>&1 | grep -o 'Installed eggNOG DB version: [^ /]*' | sed 's/.*: //')
    END_VERSIONS
    """

    stub:
    """
    touch plasmids_eggnog.emapper.annotations plasmids_eggnog.emapper.hits
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: stub
    END_VERSIONS
    """
}
