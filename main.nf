#!/usr/bin/env nextflow
/*
 * accessory_plasmids: detect accessory plasmids in Sinorhizobium meliloti assemblies.
 *
 * Input : samplesheet CSV with columns `sample,assembly` (one genome assembly per row,
 *         e.g. produced by nf-core/bacass 2.6.0).
 */

include { SEQKIT_FILTER_CONTIGS } from './modules/local/seqkit_filter_contigs'
include { GENOMAD_ENDTOEND      } from './modules/local/genomad_endtoend'

workflow {

    if (!params.input) {
        error "Please provide a samplesheet with --input (columns: sample,assembly)"
    }
    if (!params.genomad_db) {
        error "Please provide the geNomad database directory with --genomad_db"
    }

    // ---------------------------------------------------------------------
    // Read and validate the samplesheet
    // ---------------------------------------------------------------------
    ch_assemblies = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.assembly) {
                error "Samplesheet row is missing 'sample' or 'assembly': ${row}"
            }
            if (row.sample =~ /\s/) {
                error "Sample names must not contain whitespace: '${row.sample}'"
            }
            [ [id: row.sample], file(row.assembly, checkIfExists: true) ]
        }

    // Fail early on duplicated sample names (they would collide in contig IDs and outputs)
    ch_assemblies
        .map { meta, _fasta -> meta.id }
        .collect()
        .map { ids ->
            def dups = ids.countBy { it }.findAll { it.value > 1 }.keySet()
            if (dups) error "Duplicated sample names in samplesheet: ${dups.join(', ')}"
        }

    // ---------------------------------------------------------------------
    // Step 1: keep contigs <= max_contig_len (chromosome/pSymA/pSymB are removed)
    // ---------------------------------------------------------------------
    SEQKIT_FILTER_CONTIGS (
        ch_assemblies,
        params.max_contig_len,
        params.add_sample_suffix
    )

    // Genomes with no contig under the threshold yield an empty FASTA: report and drop them
    ch_small_contigs = SEQKIT_FILTER_CONTIGS.out.fasta
        .filter { meta, fasta ->
            if (fasta.size() == 0) {
                log.warn "[${meta.id}] no contigs <= ${params.max_contig_len} bp; sample skipped in downstream steps"
                return false
            }
            return true
        }

    // ---------------------------------------------------------------------
    // Step 2: plasmid / virus / chromosome classification with geNomad
    // ---------------------------------------------------------------------
    // Database is passed as a staged input so containers can see it
    ch_genomad_db = file(params.genomad_db, checkIfExists: true)

    GENOMAD_ENDTOEND (
        ch_small_contigs,
        ch_genomad_db
    )

    // geNomad is annotation only: ch_small_contigs (all filtered contigs) feeds
    // the next step; GENOMAD_ENDTOEND.out.plasmid_summary is joined later as annotation
}
