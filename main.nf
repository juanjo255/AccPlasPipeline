#!/usr/bin/env nextflow
/*
 * accessory_plasmids: detect accessory plasmids in Sinorhizobium meliloti assemblies.
 *
 * Input : samplesheet CSV with columns `sample,assembly,assembly_info` and optional
 *         `dataset` (one Flye assembly per row; assembly_info = Flye's
 *         assembly_info.txt; dataset = group used to colour plots).
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
include { REPABC_HMM            } from './modules/local/repabc_hmm'
include { CONTIG_MARKERS        } from './modules/local/contig_markers'
include { ASSEMBLY_GRAPH        } from './modules/local/assembly_graph'
include { CLEAN_GENOMES         } from './modules/local/clean_genomes'
include { BAKTA                 } from './modules/local/bakta'
include { MOB_TYPER             } from './modules/local/mob_typer'
include { PLASMID_FEATURES      } from './modules/local/plasmid_features'
include { REPABC_OPERONS        } from './modules/local/repabc_operons'
include { REPLICON_GROUPS       } from './modules/local/replicon_groups'
include { MOBMESS_SYSTEMS as MOBMESS_SYSTEMS_FINAL } from './modules/local/mobmess_systems'
include { PLOT_PLASMID_FEATURES } from './modules/local/plot_plasmid_features'
include { BAKTA_PLASMID         } from './modules/local/bakta_plasmid'
include { MERGE_ANNOTATIONS     } from './modules/local/merge_annotations'
include { TASMANIA              } from './modules/local/tasmania'
include { DEFENSEFINDER         } from './modules/local/defensefinder'
include { ANTISMASH             } from './modules/local/antismash'
include { DBCAN                 } from './modules/local/dbcan'
include { SIGNALP               } from './modules/local/signalp'
include { ANNOTATION_SUMMARY    } from './modules/local/annotation_summary'
include { TAG_REPLICONS         } from './modules/local/tag_replicons'
include { COPY_NUMBER           } from './modules/local/copy_number'
include { KMER_COMPOSITION      } from './modules/local/kmer_composition'
include { SKANI_ANI             } from './modules/local/skani_ani'
include { EGGNOG_MAPPER         } from './modules/local/eggnog_mapper'
include { COG_ANALYSIS          } from './modules/local/cog_analysis'
include { PANAROO               } from './modules/local/panaroo'
include { GENE_SHARING          } from './modules/local/gene_sharing'
include { PANAROO as PANAROO_PLASMIDS } from './modules/local/panaroo'
include { SIGNALP4              } from './modules/local/signalp4'
include { SIGNAL_PEPTIDES       } from './modules/local/signal_peptides'
include { COINTEGRATES          } from './modules/local/cointegrates'

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
    if (!params.skip_annotation) {
        def missing = ['tasmania_hmms', 'defensefinder_models', 'antismash_db', 'dbcan_db'].findAll { !params[it] }
        if (missing) {
            error "Step 10 (functional annotation) needs: ${missing.collect { '--' + it }.join(', ')} (or use --skip_annotation)"
        }
        if (params.signalp_module && !params.signalp_model_dir) {
            error "--signalp_module requires --signalp_model_dir"
        }
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
            // Optional dataset label (e.g. sequencing batch or collection) for plots; PLSDB is 'Database'
            def dataset = row.dataset ?: 'Sample'
            if (dataset == 'Database') {
                error "Dataset name 'Database' is reserved for PLSDB plasmids (sample ${row.sample})"
            }
            // Optional reads mapped back to the assembly (BAM), for plasmid copy number
            def bam = row.bam ? file(row.bam, checkIfExists: true) : []
            // Optional whole-genome GFF3 (Bakta) of the same assembly, for gene sharing with the main replicons
            def gff = row.gff ? file(row.gff, checkIfExists: true) : []
            // Optional Flye assembly graph (assembly_graph.gfa), for the graph context of linear contigs
            def graph = row.assembly_graph ? file(row.assembly_graph, checkIfExists: true) : []
            [ [id: row.sample, dataset: dataset], file(row.assembly, checkIfExists: true), file(row.assembly_info, checkIfExists: true), bam, gff, graph ]
        }

    ch_assemblies    = ch_samples.map { meta, assembly, _info, _bam, _gff, _graph -> [ meta, assembly ] }
    ch_assembly_info = ch_samples.map { meta, _assembly, info, _bam, _gff, _graph -> [ meta, info ] }
    ch_bam           = ch_samples.filter { _m, _a, _i, bam, _g, _gr -> bam }.map { meta, _a, _i, bam, _g, _gr -> [ meta, bam ] }
    ch_genome_gff    = ch_samples.filter { _m, _a, _i, _b, gff, _gr -> gff }.map { meta, _a, _i, _b, gff, _gr -> [ meta, gff ] }
    ch_graph         = ch_samples.filter { _m, _a, _i, _b, _g, graph -> graph }.map { meta, _a, info, _b, _g, graph -> [ meta, graph, info ] }

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

    // Contig -> sample -> dataset table (plots colour plasmids by dataset)
    ch_contig_dataset = SEQKIT_FILTER_CONTIGS.out.lengths
        .flatMap { meta, tsv ->
            tsv.readLines().findAll { it.trim() }.collect { line -> "${line.split('\t')[0]}\t${meta.id}\t${meta.dataset}" }
        }
        .collectFile(name: 'contig_datasets.tsv', newLine: true, sort: true, seed: 'contig\tsample\tdataset')

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

    // RepABC genes (HMMs) on every filtered contig, and circular-contig counts per assembly
    REPABC_HMM (
        ch_small_contigs
            .join(FLYE_CIRCULARITY.out.circular)
            .join(ch_assembly_info),
        file(params.repabc_hmms, checkIfExists: true)
    )

    // Flye graph context of every contig (samples with an assembly_graph column)
    ASSEMBLY_GRAPH (
        ch_graph,
        params.add_sample_suffix
    )

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

    // ---------------------------------------------------------------------
    // Step 7: RepABC operons and replicon groups of the accessory plasmids
    // ---------------------------------------------------------------------
    REPABC_OPERONS (
        BAKTA.out.tsv,
        MOBMESS_FILTER_CLUSTERS.out.clusters,
        PLASMID_FEATURES.out.ids
    )

    REPLICON_GROUPS (
        REPABC_OPERONS.out.repa,
        BAKTA.out.faa,
        file(params.repa_references, checkIfExists: true),
        params.replicon_names ? file(params.replicon_names, checkIfExists: true) : []
    )

    // ---------------------------------------------------------------------
    // Step 8: MobMess again on the final accessory plasmids (cluster IDs change)
    // ---------------------------------------------------------------------
    MOBMESS_SYSTEMS_FINAL (
        PLASMID_FEATURES.out.fasta,
        PLASMID_FEATURES.out.circular
    )

    // ---------------------------------------------------------------------
    // Step 9: summary table and figures, coloured by samplesheet dataset
    // ---------------------------------------------------------------------
    PLOT_PLASMID_FEATURES (
        PLASMID_FEATURES.out.features,
        PLASMID_FEATURES.out.rep_proteins,
        MOB_TYPER.out.mob_type,
        REPABC_OPERONS.out.summary,
        REPLICON_GROUPS.out.groups,
        ch_contig_dataset
    )

    // ---------------------------------------------------------------------
    // Step 10: functional annotation of the accessory plasmids
    //   Bakta per plasmid, then TA (TASmania), defense (DefenseFinder),
    //   BGC (antiSMASH), CAZymes (dbCAN, optional SignalP)
    // ---------------------------------------------------------------------
    if (!params.skip_annotation) {
        ch_plasmid_fasta = PLASMID_FEATURES.out.fasta
            .splitFasta(by: 1, file: true)
            .map { f -> [ [id: f.readLines()[0].substring(1).trim().split(/\s+/)[0]], f ] }

        BAKTA_PLASMID (
            ch_plasmid_fasta,
            file(params.bakta_db, checkIfExists: true)
        )

        MERGE_ANNOTATIONS (
            BAKTA_PLASMID.out.tsv.map { _meta, f -> f }.collect(),
            BAKTA_PLASMID.out.faa.map { _meta, f -> f }.collect(),
            BAKTA_PLASMID.out.gbff.map { _meta, f -> f }.collect()
        )

        // TASer-like: own Prodigal genes on the plasmid sequences, linked to Bakta genes
        TASMANIA (
            PLASMID_FEATURES.out.fasta,
            MERGE_ANNOTATIONS.out.tsv,
            file(params.tasmania_hmms, checkIfExists: true),
            params.ta_evalue,
            params.ta_report_evalue
        )

        DEFENSEFINDER (
            BAKTA_PLASMID.out.faa,
            file(params.defensefinder_models, checkIfExists: true)
        )

        ANTISMASH (
            MERGE_ANNOTATIONS.out.gbff,
            file(params.antismash_db, checkIfExists: true)
        )

        DBCAN (
            MERGE_ANNOTATIONS.out.faa,
            file(params.dbcan_db, checkIfExists: true),
            params.dbcan_methods,
            params.cazyme_min_tools
        )

        // SignalP is licensed: only run when its environment module is given
        if (params.signalp_module) {
            SIGNALP ( DBCAN.out.faa, params.signalp_model_dir )
            ch_signalp = SIGNALP.out.predictions
        } else {
            ch_signalp = channel.value([])
        }

        ANNOTATION_SUMMARY (
            PLOT_PLASMID_FEATURES.out.summary,
            TASMANIA.out.ta,
            DEFENSEFINDER.out.systems.map { _meta, f -> f }.collect(),
            ANTISMASH.out.outdir,
            DBCAN.out.cazymes,
            ch_signalp
        )
    }

    // ---------------------------------------------------------------------
    // Step 12: pairwise ANI / alignment fraction between accessory plasmids
    // ---------------------------------------------------------------------
    SKANI_ANI ( PLASMID_FEATURES.out.fasta )

    // ---------------------------------------------------------------------
    // Step 11: replicon tags (Chromosome / pSymA / pSymB), plasmid copy number
    //          and 3-mer composition vs the main replicons
    // ---------------------------------------------------------------------
    if (params.reference_genome) {
        TAG_REPLICONS (
            ch_assemblies,
            file(params.reference_genome, checkIfExists: true),
            file(params.reference_tags, checkIfExists: true),
            params.add_sample_suffix,
            params.tag_min_ref_cov
        )

        // Copy number only for samples with a BAM
        COPY_NUMBER (
            TAG_REPLICONS.out.tags.join(ch_bam),
            PLASMID_FEATURES.out.ids,
            params.add_sample_suffix
        )
        COPY_NUMBER.out.copy_number
            .map { _meta, tsv -> tsv }
            .collectFile(name: 'copy_number.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/11_replicons")

        KMER_COMPOSITION (
            TAG_REPLICONS.out.main.map { _meta, f -> f }.collect(),
            TAG_REPLICONS.out.tags.map { _meta, f -> f }.collect(),
            PLASMID_FEATURES.out.fasta,
            PLSDB_SELECT.out.ids
        )
    }

    // ---------------------------------------------------------------------
    // Step 4 (review): markers of every filtered sample contig, and linear contigs with
    // RepABC genes dropped by the MobMess filter. Placed here because the replicon tags
    // (Step 11, with --reference_genome) are used to leave out chromosome/pSym fragments.
    // ---------------------------------------------------------------------
    ch_review_tags = params.reference_genome
        ? TAG_REPLICONS.out.tags.map { _meta, f -> f }.collect()
        : channel.value([])

    CONTIG_MARKERS (
        ch_contig_dataset,
        REPABC_HMM.out.contig_info.map { _meta, f -> f }.collect(),
        REPABC_HMM.out.genes.map { _meta, f -> f }.collect(),
        REPABC_HMM.out.operons.map { _meta, f -> f }.collect(),
        REPABC_HMM.out.assembly_stats.map { _meta, f -> f }.collect(),
        GENOMAD_ENDTOEND.out.plasmid_summary.map { _meta, f -> f }.collect(),
        FLYE_CIRCULARITY.out.coverage.map { _meta, f -> f }.collect(),
        MOBMESS_SYSTEMS.out.clusters,
        MOBMESS_FILTER_CLUSTERS.out.ids,
        ch_review_tags,
        ASSEMBLY_GRAPH.out.graph.map { _meta, f -> f }.collect().ifEmpty([]),
        params.review_genomad_score
    )

    // ---------------------------------------------------------------------
    // Step 6b: clean genomes (main replicons + retained accessory plasmids + unplaced pieces),
    // re-oriented and with NCBI source modifiers. Needs the replicon tags (--reference_genome).
    // ---------------------------------------------------------------------
    if (params.reference_genome) {
        CLEAN_GENOMES (
            ch_assemblies.join(ch_assembly_info).join(TAG_REPLICONS.out.tags),
            PLASMID_FEATURES.out.ids.first(),      // value channels: shared by every sample
            CONTIG_MARKERS.out.review.first(),
            params.curation ? file(params.curation, checkIfExists: true) : [],
            file(params.repabc_hmms, checkIfExists: true),
            params.add_sample_suffix
        )
        CLEAN_GENOMES.out.summary
            .collectFile(name: 'clean_genomes_summary.tsv', keepHeader: true, skip: 1, sort: true,
                         storeDir: "${params.outdir}/06b_clean_genomes")
    }

    // ---------------------------------------------------------------------
    // Step 13: gene content (needs Step 10 annotation)
    //   A. eggNOG-mapper COG analysis of the accessory plasmids (--eggnog_db)
    //   B. Panaroo gene sharing with the main replicons (--reference_genome + samplesheet gff)
    //   C. signal peptide families (--signalp4_module)
    // ---------------------------------------------------------------------
    if (!params.skip_annotation && params.eggnog_db) {
        EGGNOG_MAPPER (
            MERGE_ANNOTATIONS.out.faa,
            file(params.eggnog_db, checkIfExists: true)
        )
        COG_ANALYSIS (
            EGGNOG_MAPPER.out.annotations,
            PLOT_PLASMID_FEATURES.out.summary,
            params.cog_min_cluster
        )
        ch_eggnog = EGGNOG_MAPPER.out.annotations
    } else {
        ch_eggnog = channel.value([])
    }

    if (!params.skip_annotation && params.reference_genome && params.gene_sharing) {
        // Genome GFF name (Panaroo sample) -> sample, to match contigs without the _<sample> suffix
        ch_gff_samples = ch_genome_gff
            .map { meta, gff -> "${gff.baseName}\t${meta.id}" }
            .collectFile(name: 'gff_samples.tsv', newLine: true, sort: true)

        PANAROO (
            ch_genome_gff.map { _meta, gff -> gff }
                .mix(BAKTA_PLASMID.out.gff3.map { _meta, gff -> gff })
                .collect()
        )
        GENE_SHARING (
            PANAROO.out.gpa,
            PANAROO.out.gene_data,
            TAG_REPLICONS.out.tags.map { _meta, f -> f }.collect(),
            PLASMID_FEATURES.out.ids,
            ch_eggnog,
            ch_gff_samples
        )
    }

    // C. Signal peptide-bearing protein families (Fig. 9D): SignalP 4.1 on all plasmid
    //    proteins + plasmid-only Panaroo families. SignalP is licensed: environment module only.
    if (!params.skip_annotation && params.signalp4_module) {
        PANAROO_PLASMIDS (
            BAKTA_PLASMID.out.gff3.map { _meta, gff -> gff }.collect()
        )
        SIGNALP4 ( MERGE_ANNOTATIONS.out.faa )
        SIGNAL_PEPTIDES (
            SIGNALP4.out.summary,
            MERGE_ANNOTATIONS.out.faa,
            PANAROO_PLASMIDS.out.gpa,
            PLOT_PLASMID_FEATURES.out.summary,
            DBCAN.out.cazymes,
            params.signalp_top
        )
    }

    // ---------------------------------------------------------------------
    // Step 14: cointegrates (MobMess containment of Step 4 + Bakta junction genes)
    // ---------------------------------------------------------------------
    if (!params.skip_annotation) {
        COINTEGRATES (
            MOBMESS_SYSTEMS.out.clusters,
            MOBMESS_SYSTEMS.out.systems,
            MOBMESS_SYSTEMS.out.ani_blocks,
            PLASMID_FEATURES.out.ids,
            MERGE_ANNOTATIONS.out.tsv
        )
    }
}
