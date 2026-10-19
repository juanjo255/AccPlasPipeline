/*
 * Signal peptide-bearing protein families on the accessory plasmids (Fig. 9D):
 * SignalP 4.1 hits linked to the plasmid-only Panaroo families (bin/signal_peptides.py).
 */
process SIGNAL_PEPTIDES {
    tag "all"
    label 'process_single'

    conda "${moduleDir}/environment.yml"

    input:
    path signalp
    path proteins
    path gpa
    path plasmid_summary
    path cazymes
    val top

    output:
    path "signal_peptide_*.tsv"      , emit: tables
    path "signal_peptide_summary.tsv", emit: summary
    path "*.{png,svg}"               , emit: figures
    path "versions.yml"              , emit: versions

    script:
    def caz = cazymes ? "--cazymes ${cazymes}" : ''
    """
    signal_peptides.py \\
        --signalp ${signalp} \\
        --proteins ${proteins} \\
        --gpa ${gpa} \\
        --plasmids ${plasmid_summary} \\
        ${caz} \\
        --top ${top}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch signal_peptide_proteins.tsv signal_peptide_families.tsv signal_peptide_summary.tsv \\
        signal_peptide_families.png signal_peptide_families.svg
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
