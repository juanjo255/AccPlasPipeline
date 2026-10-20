/*
 * Clean genome of one strain (bin/clean_genome.py): chromosome, pSymA, pSymB, retained accessory
 * plasmids and unplaced replicon pieces / review contigs, with NCBI source modifiers. The circular
 * chromosome is re-oriented to dnaA (dnaapler), circular plasmids to the RepA of their RepABC
 * operon (RepABC HMMs). Everything left out is written to <sample>.removed.fasta with the reason.
 */
process CLEAN_GENOMES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(assembly), path(assembly_info), path(tags)
    path plasmid_ids
    path review
    path curation
    path hmm
    val add_sample_suffix

    output:
    tuple val(meta), path("${meta.id}.clean.fasta")  , emit: fasta
    tuple val(meta), path("${meta.id}.clean.tsv")    , emit: table
    tuple val(meta), path("${meta.id}.removed.fasta"), emit: removed
    tuple val(meta), path("${meta.id}.rotation.tsv") , emit: rotation
    path "${meta.id}.summary.tsv"                    , emit: summary
    path "versions.yml"                              , emit: versions

    script:
    def suffix   = add_sample_suffix ? '--add-suffix' : ''
    def rev_arg  = review   ? "--review ${review}"     : ''
    def cur_arg  = curation ? "--curation ${curation}" : ''
    """
    clean_genome.py select \\
        --assembly ${assembly} \\
        --assembly-info ${assembly_info} \\
        --tags ${tags} \\
        --plasmids ${plasmid_ids} \\
        ${rev_arg} ${cur_arg} \\
        --sample ${meta.id} ${suffix}

    # Circular chromosome -> dnaapler (dnaA); linear chromosomes keep their orientation
    chrom_arg=''
    chrom=\$(awk -F'\\t' '\$4=="chromosome" && \$6=="circular" {print \$1}' ${meta.id}.clean.tsv | head -1)
    if [ -n "\$chrom" ]; then
        awk -v c=">\$chrom" '/^>/{p=(\$1==c)} p' ${meta.id}.clean.plain.fasta > chromosome.fasta
        dnaapler chromosome -i chromosome.fasta -o dnaapler -p ${meta.id} -t ${task.cpus} -f
        chrom_arg="--chromosome dnaapler/${meta.id}_reoriented.fasta"
    fi

    # RepABC on the clean sequences -> RepA of each circular plasmid
    awk -F'\\t' 'NR>1 && \$4!="removed" {print \$1"\\t"(\$6=="circular"?1:0)}' ${meta.id}.clean.tsv > circular.tsv
    repabc_hmm.py --fasta ${meta.id}.clean.plain.fasta --hmm ${hmm} --circular circular.tsv \\
        --threads ${task.cpus} --prefix repabc
    repabc_operons.py --bakta repabc.repabc.tsv --clusters repabc.contig_info.tsv --prefix repabc

    clean_genome.py rotate \\
        --fasta ${meta.id}.clean.plain.fasta \\
        --table ${meta.id}.clean.tsv \\
        --repabc repabc.repabc.tsv \\
        --operons repabc.operons.tsv \\
        \$chrom_arg \\
        --sample ${meta.id}

    clean_genome.py finalize \\
        --fasta ${meta.id}.reoriented.fasta \\
        --table ${meta.id}.clean.tsv \\
        --sample ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$(dnaapler --version | sed 's/.*version //')
        pyhmmer: \$(python3 -c 'import pyhmmer; print(pyhmmer.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.clean.fasta ${meta.id}.clean.tsv ${meta.id}.removed.fasta ${meta.id}.rotation.tsv
    printf 'sample\\tchromosome_bp\\tpSymA_bp\\tpSymB_bp\\taccessory_plasmids\\taccessory_circular\\tunplaced\\tremoved\\ttotal_bp\\tnotes\\n${meta.id}\\t-\\t-\\t-\\t0\\t0\\t0\\t0\\t0\\tstub\\n' > ${meta.id}.summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: stub
    END_VERSIONS
    """
}
