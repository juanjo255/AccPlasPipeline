#!/usr/bin/env nextflow
/*
 * accessory_plasmids: detect accessory plasmids in Sinorhizobium meliloti assemblies.
 *
 * Input : samplesheet CSV with columns `sample,assembly,assembly_info` (one Flye
 *         assembly per row; assembly_info = Flye's assembly_info.txt).
 */

include { SEQKIT_FILTER_CONTIGS } from './modules/local/seqkit_filter_contigs'
include { GENOMAD_ENDTOEND      } from './modules/local/genomad_endtoend'
include { MASH_DIST_PLSDB       } from './modules/local/mash_dist_plsdb'
include { PLSDB_SELECT          } from './modules/local/plsdb_select'
include { PLSDB_EXTRACT         } from './modules/local/plsdb_extract'
include { FLYE_CIRCULARITY      } from './modules/local/flye_circularity'
include { MOBMESS_SYSTEMS       } from './modules/local/mobmess_systems'
include { MOBMESS_FILTER_CLUSTERS } from './modules/local/mobmess_filter_clusters'
include { MOBMESS_NETWORK       } from './modules/local/mobmess_network'
include { BAKTA                 } from './modules/local/bakta'
include { MOB_TYPER             } from './modules/local/mob_typer'
include { PLASMID_FEATURES      } from './modules/local/plasmid_features'

workflow {

    if (!params.input) {
        error "Please provide a samplesheet with --input (columns: sample,assembly,assembly_info)"
    }
    if (!params.genomad_db) {
        error "Please provide the geNomad database directory with --genomad_db"
    }
    if (!params.plsdb_fasta || !params.plsdb_nuccore || !params.plsdb_typing) {
        error "Please provide the PLSDB files --plsdb_fasta, --plsdb_nuccore (nuccore.csv) and --plsdb_typing (typing.csv)"
    }
    if (!params.bakta_db || !params.mobsuite_db) {
        error "Please provide the Bakta database (--bakta_db) and MOB-suite database (--mobsuite_db)"
    }
    if (params.plsdb_mash && !params.plsdb_sketch) {
        error "--plsdb_mash requires the PLSDB Mash sketch (--plsdb_sketch)"
    }

    // ---------------------------------------------------------------------
    // Read and validate the samplesheet
    // ---------------------------------------------------------------------
    ch_samples = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.assembly || !row.assembly_info) {
                error "Samplesheet row is missing 'sample', 'assembly' or 'assembly_info': ${row}"
            }
            if (row.sample =~ /\s/) {
                error "Sample names must not contain whitespace: '${row.sample}'"
            }
            [ [id: row.sample], file(row.assembly, checkIfExists: true), file(row.assembly_info, checkIfExists: true) ]
        }

    ch_assemblies    = ch_samples.map { meta, assembly, _info -> [ meta, assembly ] }
    ch_assembly_info = ch_samples.map { meta, _assembly, info -> [ meta, info ] }

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

    // ---------------------------------------------------------------------
    // Step 3: related PLSDB plasmids for de-replication
    //   metadata rule (always) + mash hits to the filtered contigs (optional)
    // ---------------------------------------------------------------------
    if (params.plsdb_mash) {
        MASH_DIST_PLSDB (
            ch_small_contigs,
            file(params.plsdb_sketch, checkIfExists: true),
            params.mash_max_dist,
            params.mash_max_pvalue
        )
        ch_mash_hits = MASH_DIST_PLSDB.out.dist.map { _meta, tsv -> tsv }.collect()
    } else {
        ch_mash_hits = channel.value([])
    }

    PLSDB_SELECT (
        file(params.plsdb_nuccore, checkIfExists: true),
        file(params.plsdb_typing, checkIfExists: true),
        params.plsdb_keyword,
        ch_mash_hits
    )

    PLSDB_EXTRACT (
        file(params.plsdb_fasta, checkIfExists: true),
        PLSDB_SELECT.out.ids
    )

    // ---------------------------------------------------------------------
    // Step 4: de-replication with MobMess
    // ---------------------------------------------------------------------
    // Circularity of the sample contigs from Flye's assembly_info.txt
    FLYE_CIRCULARITY (
        ch_small_contigs.join(ch_assembly_info),
        params.add_sample_suffix
    )

    // Pool sample contigs + related PLSDB plasmids (and their circularity) into one run
    ch_mobmess_fasta = ch_small_contigs
        .map { _meta, fasta -> fasta }
        .mix(PLSDB_EXTRACT.out.fasta)
        .collect()
    ch_mobmess_circular = FLYE_CIRCULARITY.out.circular
        .map { _meta, tsv -> tsv }
        .mix(PLSDB_SELECT.out.circular)
        .collect()

    MOBMESS_SYSTEMS (
        ch_mobmess_fasta,
        ch_mobmess_circular
    )

    // Keep clusters that look like real plasmids (see README, Step 4)
    MOBMESS_FILTER_CLUSTERS (
        MOBMESS_SYSTEMS.out.clusters,
        MOBMESS_SYSTEMS.out.fasta,
        PLSDB_SELECT.out.ids,
        params.mobmess_max_len
    )

    // ---------------------------------------------------------------------
    // Step 5: MobMess containment network plot (all clusters, before filtering)
    // ---------------------------------------------------------------------
    MOBMESS_NETWORK (
        MOBMESS_SYSTEMS.out.contig_graph,
        MOBMESS_SYSTEMS.out.fasta,
        PLSDB_SELECT.out.ids,
        MOBMESS_FILTER_CLUSTERS.out.ids,
        GENOMAD_ENDTOEND.out.plasmid_summary.map { _meta, tsv -> tsv }.collect()
    )

    // ---------------------------------------------------------------------
    // Step 6: plasmid markers on sequences that passed the MobMess filter
    //   Bakta (rep/partition proteins), geNomad (MOB relaxase), MOB-typer (report)
    // ---------------------------------------------------------------------
    BAKTA (
        MOBMESS_FILTER_CLUSTERS.out.fasta,
        file(params.bakta_db, checkIfExists: true)
    )

    MOB_TYPER (
        MOBMESS_FILTER_CLUSTERS.out.fasta,
        file(params.mobsuite_db, checkIfExists: true)
    )

    PLASMID_FEATURES (
        MOBMESS_FILTER_CLUSTERS.out.fasta,
        MOBMESS_FILTER_CLUSTERS.out.clusters,
        BAKTA.out.tsv,
        MOB_TYPER.out.mob_type,
        GENOMAD_ENDTOEND.out.plasmid_summary.map { _meta, tsv -> tsv }.collect(),
        params.feature_exclude ? file(params.feature_exclude, checkIfExists: true) : [],
        FLYE_CIRCULARITY.out.coverage.map { _meta, tsv -> tsv }.collect(),
        params.flag_max_len,
        params.flag_min_rel_cov
    )

    // PLASMID_FEATURES.out.fasta (final accessory plasmids) feeds the next step of the pipeline
}
