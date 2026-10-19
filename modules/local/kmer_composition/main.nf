/*
 * Trinucleotide composition of accessory plasmids vs chromosome / pSymA / pSymB:
 * Kmertools 3-mer frequencies, Bray-Curtis distances, UMAP, length and GC plots
 * (bin/plasmid_kmer_composition.py).
 */
process KMER_COMPOSITION {
    tag "all"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    path main_replicons, stageAs: 'main/*'
    path tags          , stageAs: 'tags/*'
    path plasmids
    path exclude

    output:
    path "kmer_3mer_frequencies.tsv", emit: frequencies
    path "kmer_umap.tsv"            , emit: umap
    path "kmer_stats.tsv"           , emit: stats
    path "*.{png,svg}"              , emit: figures
    path "versions.yml"             , emit: versions

    script:
    """
    plasmid_kmer_composition.py prepare \\
        --main main/* \\
        --tags tags/* \\
        --plasmids ${plasmids} \\
        --exclude ${exclude}

    kmertools comp oligo -k 3 -i kmer_input.fasta -o kmer_3mer.txt

    plasmid_kmer_composition.py analyse --kmers kmer_3mer.txt --table kmer_input.tsv
    rm -f kmer_input.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmertools: \$(kmertools --version | sed 's/^kmertools //')
        umap-learn: \$(python3 -c 'import importlib.metadata as m; print(m.version("umap-learn"))')
    END_VERSIONS
    """

    stub:
    """
    touch kmer_3mer_frequencies.tsv kmer_umap.tsv kmer_stats.tsv kmer_umap.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmertools: stub
    END_VERSIONS
    """
}
