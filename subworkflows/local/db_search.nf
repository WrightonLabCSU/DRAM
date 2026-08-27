//
// Subworkflow with functionality specific to the BortonWrightonLabs/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GENE_LOCS                                     } from "../../modules/local/annotate/gene_locs.nf"

include { COMBINE_ANNOTATIONS                           } from "../../modules/local/annotate/combine_annotations.nf"

include { MMSEQS_INDEX                                  } from "../../modules/local/annotate/mmseqs_index.nf"

// NextFlow only process with the same name in the same workflow, so either alias it or include it a different workflow
include { MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS         } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL          } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER         } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL         } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD        } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG           } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF         } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM           } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CARD           } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_TCDB           } from "../../modules/local/annotate/mmseqs_search.nf"

include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_MEROPS } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_VIRAL  } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_CAMPER } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_METHYL } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_CANTHYD } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_KEGG   } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_UNIREF } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_CARD   } from "../../modules/local/annotate/mmseqs_search_gpu.nf"
include { MMSEQS_SEARCH_GPU as MMSEQS_SEARCH_GPU_TCDB   } from "../../modules/local/annotate/mmseqs_search_gpu.nf"

include { ADD_SQL_DESCRIPTIONS as SQL_UNIREF            } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_VIRAL             } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_MEROPS            } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_KEGG              } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_PFAM              } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_DBCAN             } from "../../modules/local/annotate/add_sql_descriptions.nf"

include { HMM_SEARCH as HMM_SEARCH_KOFAM                } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_DRAM_DB              } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_VOG                  } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CAMPER               } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CANTHYD              } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_SULFUR               } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_FEGENIE              } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_METALS               } from "../../modules/local/annotate/hmmsearch.nf"

include { ANTISMASH_ANTISMASH                           } from '../../modules/nf-core/antismash/antismash/main'
include { RGI_MAIN                                      } from '../../modules/nf-core/rgi/main/main'
include { RUNDBCAN_EASYSUBSTRATE                        } from '../../modules/nf-core/rundbcan/easysubstrate/main'

include {checkDBVersion                                 } from '../../subworkflows/local/utils_pipeline_setup.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO DB_SEARCH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DB_SEARCH {
    take:
    ch_gene_locs  // channel: path(gene_locs_tsv) ]
    ch_called_proteins  // channel: [ val(input_fasta name), path(called_proteins file) ]
    ch_filtered_fasta  // channel: [ val(input_fasta name), path(filtered_fasta file) ]
    ch_gene_gff  // channel: [ val(input_fasta name), path(gene_gff file) ]
    ch_called_genes  // channel: [ val(input_fasta name), path(called_genes file) ]
    default_sheet // Path to dummy sheet
    use_kegg
    use_kofam
    use_dbcan
    use_camper
    use_fegenie
    use_methyl
    use_canthyd
    use_sulfur
    use_pfam
    use_merops
    use_uniref
    use_metals
    use_antismash
    use_rgi
    use_card
    use_tcdb
    use_dram_db
    use_vog

    main:

    DB_CHANNEL_SETUP(
        use_kegg,
        use_kofam,
        use_dbcan,
        use_camper,
        use_fegenie,
        use_methyl,
        use_canthyd,
        use_sulfur,
        use_pfam,
        use_merops,
        use_uniref,
        use_metals,
        use_antismash,
        use_rgi,
        use_card,
        use_tcdb,
        use_dram_db,
        use_vog

    )

    ch_sql_descriptions_db = file(params.sql_descriptions_db)
    ch_kofam_list = file(params.kofam_list)
    ch_vog_list = file(params.vog_list)
    ch_camper_hmm_list = file(params.camper_hmm_list)
    ch_canthyd_hmm_list = file(params.cant_hyd_hmm_list)
    ch_dram_db_hmm_list = file(params.dram_db_list)


    kegg_name = "kegg"
    dbcan_name = "dbcan"
    kofam_name = "kofam"
    merops_name = "merops"
    viral_name = "viral"
    camper_name = "camper"
    canthyd_name = "cant_hyd"
    fegenie_name = "fegenie"
    sulfur_name = "sulfur"
    methyl_name = "methyl"
    uniref_name = "uniref"
    pfam_name = "pfam"
    vogdb_name = "vogdb"
    metals_name = "metals"
    card_name = "card"
    tcdb_name = "tcdb"
    dram_db_name = "dram_db"

    def formattedOutputChannels = channel.of()
    def dbcanOutputChannels = channel.of()
    use_mmseqs_gpu = workflow.profile.contains('gpu')
    mmseqs_gpu_excluded_dbs = (params.mmseqs_gpu_exclude_dbs?.tokenize(',')?.collect { db -> db.trim().toLowerCase() } ?: []).findAll { db -> db && db != 'none' }
    use_gpu_kegg = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('kegg')
    use_gpu_camper = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('camper')
    use_gpu_methyl = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('methyl')
    use_gpu_canthyd = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('canthyd')
    use_gpu_merops = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('merops')
    use_gpu_uniref = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('uniref')
    use_gpu_card = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('card')
    use_gpu_tcdb = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('tcdb')
    use_gpu_viral = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('viral')
    // Here we will create mmseqs2 index files for each of the inputs if we are going to do a mmseqs2 database
    // We use .val because we need to unwrap the workflow output.
    // if the .out was from a process, this could block as it waited for DB_CHANNEL_SETUP, so use with caution
    if (DB_CHANNEL_SETUP.out.index_mmseqs.val) {
        // Use MMSEQS2 to index each called genes protein file
        MMSEQS_INDEX( ch_called_proteins )
        ch_mmseqs_query = MMSEQS_INDEX.out.mmseqs_index_out
    }

    // KEGG annotation
    if (use_kegg) {
        ch_combined_query_locs_kegg = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_kegg) {
            MMSEQS_SEARCH_GPU_KEGG( ch_combined_query_locs_kegg, DB_CHANNEL_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, kegg_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_GPU_KEGG.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_KEGG( ch_combined_query_locs_kegg, DB_CHANNEL_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, kegg_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_KEGG.out.mmseqs_search_formatted_out
        }

        SQL_KEGG(ch_mmseqs_unformatted, kegg_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_KEGG.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // KOFAM annotation
    if (use_kofam) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_KOFAM (
            ch_combined_proteins_locs,
            params.kofam_e_value,
            DB_CHANNEL_SETUP.out.ch_kofam_db,
            ch_kofam_list,
            true,
            kofam_name
            )
        ch_hmm_formatted = HMM_SEARCH_KOFAM.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // PFAM annotation
    if (use_pfam) {
        ch_combined_query_locs_pfam = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_PFAM( ch_combined_query_locs_pfam, DB_CHANNEL_SETUP.out.ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, pfam_name )
        ch_mmseqs_unformatted = MMSEQS_SEARCH_PFAM.out.mmseqs_search_formatted_out

        SQL_PFAM(ch_mmseqs_unformatted, pfam_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_PFAM.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }

    // dbCAN3 annotation
    if  (use_dbcan) {

        ch_dbcan_inputs = ch_called_proteins
            .join(ch_gene_gff, by: [0])
            .multiMap { name, called_proteins, gene_gff ->
                proteins:
                    tuple([id: name], called_proteins)

                gff:
                    tuple([id: name], gene_gff, "prodigal")
            }
        RUNDBCAN_EASYSUBSTRATE(
            ch_dbcan_inputs.proteins,
            ch_dbcan_inputs.gff,
            DB_CHANNEL_SETUP.out.ch_dbcan_db
        )
        dbcanOutputChannels = dbcanOutputChannels.mix(RUNDBCAN_EASYSUBSTRATE.out.dbcanhmm_results)
        dbcanOutputChannels = dbcanOutputChannels.mix(RUNDBCAN_EASYSUBSTRATE.out.dbcansub_results)
    }
    // CAMPER annotation
    if (use_camper) {
        // HMM
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_CAMPER (
            ch_combined_proteins_locs,
            params.camper_e_value,
            DB_CHANNEL_SETUP.out.ch_camper_hmm_db,
            ch_camper_hmm_list,
            false,
            camper_name
        )
        ch_hmm_formatted = HMM_SEARCH_CAMPER.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)

        // MMseqs
        ch_combined_query_locs_camper = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_camper) {
            MMSEQS_SEARCH_GPU_CAMPER( ch_combined_query_locs_camper, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_list, camper_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_GPU_CAMPER.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_CAMPER( ch_combined_query_locs_camper, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_list, camper_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_CAMPER.out.mmseqs_search_formatted_out
        }
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // FeGenie annotation
    if (use_fegenie) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_FEGENIE (
            ch_combined_proteins_locs,
            params.fegenie_e_value,
            DB_CHANNEL_SETUP.out.ch_fegenie_db,
            default_sheet,
            false,
            fegenie_name
            )
        ch_hmm_formatted = HMM_SEARCH_FEGENIE.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // Methyl annotation
    if (use_methyl) {
        ch_combined_query_locs_methyl = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_methyl) {
            MMSEQS_SEARCH_GPU_METHYL( ch_combined_query_locs_methyl, DB_CHANNEL_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, methyl_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_GPU_METHYL.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_METHYL( ch_combined_query_locs_methyl, DB_CHANNEL_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, methyl_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_METHYL.out.mmseqs_search_formatted_out
        }
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // CANT-HYD annotation
    if (use_canthyd) {
        // MMseqs
        ch_combined_query_locs_canthyd = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_canthyd) {
            MMSEQS_SEARCH_GPU_CANTHYD( ch_combined_query_locs_canthyd, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_GPU_CANTHYD.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_CANTHYD( ch_combined_query_locs_canthyd, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_CANTHYD.out.mmseqs_search_formatted_out
        }
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)

        //HMM
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_CANTHYD (
            ch_combined_proteins_locs,
            params.canthyd_e_value,
            DB_CHANNEL_SETUP.out.ch_canthyd_hmm_db,
            ch_canthyd_hmm_list,
            false,
            canthyd_name
            )
        ch_hmm_formatted = HMM_SEARCH_CANTHYD.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // Sulfur annotation
    if (use_sulfur) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_SULFUR (
            ch_combined_proteins_locs,
            params.sulfur_e_value,
            DB_CHANNEL_SETUP.out.ch_sulfur_db,
            default_sheet,
            false,
            sulfur_name
            )
        ch_hmm_formatted = HMM_SEARCH_SULFUR.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // MEROPS annotation
    if (use_merops) {
        ch_combined_query_locs_merops = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_merops) {
            MMSEQS_SEARCH_GPU_MEROPS( ch_combined_query_locs_merops, DB_CHANNEL_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, merops_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_GPU_MEROPS.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_MEROPS( ch_combined_query_locs_merops, DB_CHANNEL_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, merops_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_MEROPS.out.mmseqs_search_formatted_out
        }

        SQL_MEROPS(ch_mmseqs_unformatted, merops_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_MEROPS.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // Uniref annotation
    if (use_uniref) {
        ch_combined_query_locs_uniref = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_uniref) {
            MMSEQS_SEARCH_GPU_UNIREF( ch_combined_query_locs_uniref, DB_CHANNEL_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, uniref_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_GPU_UNIREF.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_UNIREF( ch_combined_query_locs_uniref, DB_CHANNEL_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, uniref_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_UNIREF.out.mmseqs_search_formatted_out
        }

        SQL_UNIREF(ch_mmseqs_unformatted, uniref_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_UNIREF.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // Metals annotation
    if (use_metals) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_METALS (
            ch_combined_proteins_locs,
            params.metals_e_value,
            DB_CHANNEL_SETUP.out.ch_metals_db,
            default_sheet,
            false,
            metals_name
            )
        ch_hmm_formatted = HMM_SEARCH_METALS.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // antiSMASH
    if (use_antismash) {
        ch_filtered_fasta.ifEmpty{ log.warn("Antismash requires raw fasta files, skipping antismash") }
        ch_antismash_inputs = ch_filtered_fasta
            .join(ch_gene_gff, by: [0])
            .multiMap { name, filtered_fasta, gene_gff ->
                fasta:
                    tuple([id: name], filtered_fasta)

                gff:
                    gene_gff
            }
        ANTISMASH_ANTISMASH(
            ch_antismash_inputs.fasta,
            DB_CHANNEL_SETUP.out.ch_antismash_db,
            ch_antismash_inputs.gff
        )
    }
    // RGI with CARD
    if (use_rgi) {
        RGI_MAIN(
            ch_called_genes.map { name, called_genes -> tuple([id: name], called_genes)},
            DB_CHANNEL_SETUP.out.ch_card_db,
            []
        )
    }
    // CARD annotation
    if (use_card) {
        ch_combined_query_locs_card = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_card) {
            MMSEQS_SEARCH_GPU_CARD( ch_combined_query_locs_card, DB_CHANNEL_SETUP.out.ch_card_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, card_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_GPU_CARD.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_CARD( ch_combined_query_locs_card, DB_CHANNEL_SETUP.out.ch_card_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, card_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_CARD.out.mmseqs_search_formatted_out
        }
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // TCDB annotation
    if (use_tcdb) {
        ch_combined_query_locs_tcdb = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_tcdb) {
            MMSEQS_SEARCH_GPU_TCDB( ch_combined_query_locs_tcdb, DB_CHANNEL_SETUP.out.ch_tcdb_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, tcdb_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_GPU_TCDB.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_TCDB( ch_combined_query_locs_tcdb, DB_CHANNEL_SETUP.out.ch_tcdb_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, tcdb_name )
            ch_mmseqs_formatted = MMSEQS_SEARCH_TCDB.out.mmseqs_search_formatted_out
        }
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // DRAM DB annotation
     if (use_dram_db) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_DRAM_DB (
            ch_combined_proteins_locs,
            "",  // No e value, skip e value flag
            DB_CHANNEL_SETUP.out.ch_dram_db,
            ch_dram_db_hmm_list,
            false,
            dram_db_name
            )
        ch_hmm_formatted = HMM_SEARCH_DRAM_DB.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // VOGdb annotation
    if (use_vog) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_VOG (
            ch_combined_proteins_locs,
            params.vog_e_value,
            DB_CHANNEL_SETUP.out.ch_vogdb_db,
            default_sheet,
            false,
            vogdb_name
            )
        ch_hmm_formatted = HMM_SEARCH_VOG.out.formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // Viral annotation
    if (params.use_viral) {
        ch_combined_query_locs_viral = ch_mmseqs_query.join(ch_gene_locs)
        if (use_gpu_viral) {
            MMSEQS_SEARCH_GPU_VIRAL( ch_combined_query_locs_viral, DB_CHANNEL_SETUP.out.ch_viral_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, viral_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_GPU_VIRAL.out.mmseqs_search_formatted_out
        } else {
            MMSEQS_SEARCH_VIRAL( ch_combined_query_locs_viral, DB_CHANNEL_SETUP.out.ch_viral_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, viral_name )
            ch_mmseqs_unformatted = MMSEQS_SEARCH_VIRAL.out.mmseqs_search_formatted_out
        }

        SQL_VIRAL(ch_mmseqs_unformatted, viral_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_VIRAL.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }

    fastas = formattedOutputChannels.map { it[1] }.toList()
    genes = ch_called_proteins.map { it[1] }.toList()
    dbcan_output = dbcanOutputChannels.map { it[1] }.toList()
    COMBINE_ANNOTATIONS(
        fastas,
        genes,
        dbcan_output
    )
    ch_combined_annotations = COMBINE_ANNOTATIONS.out.combined_annotations_out


    emit:
    ch_combined_annotations  // channel: [ path(combined_annotations_out) ]

}

workflow DB_CHANNEL_SETUP {
    take:
    use_kegg
    use_kofam
    use_dbcan
    use_camper
    use_fegenie
    use_methyl
    use_canthyd
    use_sulfur
    use_pfam
    use_merops
    use_uniref
    use_metals
    use_antismash
    use_rgi
    use_card
    use_tcdb
    use_dram_db
    use_vog


    main:


    use_mmseqs_gpu = workflow.profile.contains('gpu')
    mmseqs_gpu_excluded_dbs = (params.mmseqs_gpu_exclude_dbs?.tokenize(',')?.collect { db -> db.trim().toLowerCase() } ?: []).findAll { db -> db && db != 'none' }
    use_gpu_kegg = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('kegg')
    use_gpu_camper = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('camper')
    use_gpu_methyl = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('methyl')
    use_gpu_canthyd = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('canthyd')
    use_gpu_merops = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('merops')
    use_gpu_uniref = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('uniref')
    use_gpu_card = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('card')
    use_gpu_tcdb = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('tcdb')
    use_gpu_viral = use_mmseqs_gpu && !mmseqs_gpu_excluded_dbs.contains('viral')
    index_mmseqs = false
    ch_kegg_db = channel.empty()
    ch_kofam_db = channel.empty()
    ch_dbcan_db = channel.empty()
    ch_camper_hmm_db = channel.empty()
    ch_camper_mmseqs_db = channel.empty()
    ch_camper_mmseqs_list = channel.empty()
    ch_merops_db = channel.empty()
    ch_pfam_mmseqs_db = channel.empty()
    ch_heme_db = channel.empty()
    ch_sulfur_db = channel.empty()
    ch_uniref_db = channel.empty()
    ch_metals_db = channel.empty()
    ch_antismash_db = channel.empty()
    ch_card_db = channel.empty()
    ch_card_mmseqs_db = channel.empty()
    ch_tcdb_db = channel.empty()
    ch_dram_db = channel.empty()
    ch_methyl_db = channel.empty()
    ch_fegenie_db = channel.empty()
    ch_canthyd_hmm_db = channel.empty()
    ch_canthyd_mmseqs_db = channel.empty()
    ch_canthyd_mmseqs_list = channel.empty()
    ch_vogdb_db = channel.empty()
    ch_viral_db = channel.empty()

    if (use_kegg) {
        kegg_db_path = use_gpu_kegg ? params.kegg_gpu_db : params.kegg_db
        ch_kegg_db = file(kegg_db_path).exists() ? file(kegg_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. KEGG ${use_gpu_kegg ? 'gpu' : 'cpu'} database file not found at ${kegg_db_path}")
        index_mmseqs = true
    }

    if (use_kofam) {
        ch_kofam_db = file(params.kofam_db).exists() ? file(params.kofam_db) : error("Error: If using --annotate, you must supply prebuilt databases. KOFAM database file not found at ${params.kofam_db}")
    }

    if (use_dbcan) {
        ch_dbcan_db = file(params.dbcan_db).exists() ? file(params.dbcan_db) : error("Error: If using --annotate, you must supply prebuilt databases. DBCAN database file not found at ${params.dbcan_db}")
        checkDBVersion(params.dbcan_version_file, params.dbcan_version, "dbcan")
    }

    if (use_camper) {
        ch_camper_hmm_db = file(params.camper_hmm_db).exists() ? file(params.camper_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER HMM database file not found at ${params.camper_hmm_db}")
        camper_mmseqs_db_path = use_gpu_camper ? params.camper_mmseqs_gpu_db : params.camper_mmseqs_db
        ch_camper_mmseqs_db = file(camper_mmseqs_db_path).exists() ? file(camper_mmseqs_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER MMseqs2 ${use_gpu_camper ? 'gpu' : 'cpu'} database file not found at ${camper_mmseqs_db_path}")
        index_mmseqs = true
        ch_camper_mmseqs_list = file(params.camper_mmseqs_list)
    }

    if (use_merops) {
        merops_db_path = use_gpu_merops ? params.merops_gpu_db : params.merops_db
        ch_merops_db = file(merops_db_path).exists() ? file(merops_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. MEROPS ${use_gpu_merops ? 'gpu' : 'cpu'} database file not found at ${merops_db_path}")
        index_mmseqs = true
    }

    if (use_pfam) {
        ch_pfam_mmseqs_db = file(params.pfam_mmseq_db).exists() ? file(params.pfam_mmseq_db) : error("Error: If using --annotate, you must supply prebuilt databases. PFAM database file not found at ${params.pfam_mmseq_db}")
        index_mmseqs = true
    }

    // if (use_heme) {
    //     ch_heme_db = file(params.heme_db).exists() ? file(params.heme_db) : error("Error: If using --annotate, you must supply prebuilt databases. HEME database file not found at ${params.heme_db}")
    // }

    if (use_sulfur) {
        ch_sulfur_db = file(params.sulfur_db).exists() ? file(params.sulfur_db) : error("Error: If using --annotate, you must supply prebuilt databases. SULURR database file not found at ${params.sulfur_db}")
    }

    if (use_uniref) {
        uniref_db_path = use_gpu_uniref ? params.uniref_gpu_db : params.uniref_db
        ch_uniref_db = file(uniref_db_path).exists() ? file(uniref_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. UNIREF ${use_gpu_uniref ? 'gpu' : 'cpu'} database file not found at ${uniref_db_path}")
        index_mmseqs = true
    }

    if (use_metals) {
        ch_metals_db = file(params.metals_db).exists() ? file(params.metals_db) : error("Error: If using --annotate, you must supply prebuilt databases. METALS database file not found at ${params.metals_db}")
    }

    if (use_antismash) {
        ch_antismash_db = file(params.antismash_db).exists() ? file(params.antismash_db) : error("Error: If using --annotate, you must supply prebuilt databases. antismash database file not found at ${params.antismash_db}")
    }

    if (use_rgi) {
        ch_card_db = file(params.card_db).exists() ? file(params.card_db) : error("Error: If using --annotate, you must supply prebuilt databases. rgi database file not found at ${params.card_db}")
    }

    if (use_card) {
        card_mmseqs_db_path = use_gpu_card ? params.card_gpu_db : params.card_db
        ch_card_mmseqs_db = file(card_mmseqs_db_path).exists() ? file(card_mmseqs_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. CARD MMseqs2 ${use_gpu_card ? 'gpu' : 'cpu'} database file not found at ${card_mmseqs_db_path}")
        index_mmseqs = true
    }

    if (use_tcdb) {
        tcdb_db_path = use_gpu_tcdb ? params.tcdb_gpu_db : params.tcdb_db
        ch_tcdb_db = file(tcdb_db_path).exists() ? file(tcdb_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. TCDB ${use_gpu_tcdb ? 'gpu' : 'cpu'} database file not found at ${tcdb_db_path}")
        index_mmseqs = true
    }

    if (use_dram_db) {
        if (!file(params.dram_db).exists()) {
            error("Error: If using --annotate, you must supply prebuilt databases. dram database file not found at ${params.dram_db}")
        }
        // ch_dram_db = [file("${params.dram_db}/dram_db.hmm")]
        // ch_dram_db = [file(params.dram_db)]
        ch_dram_db = file(params.dram_db)
        // ch_dram_db = file(params.dram_db).exists() ? file(params.dram_db) : error("Error: If using --annotate, you must supply prebuilt databases. dram database file not found at ${params.dram_db}")
    }

    if (use_methyl) {
        methyl_db_path = use_gpu_methyl ? params.methyl_gpu_db : params.methyl_db
        ch_methyl_db = file(methyl_db_path).exists() ? file(methyl_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. METHYL ${use_gpu_methyl ? 'gpu' : 'cpu'} database file not found at ${methyl_db_path}")
        index_mmseqs = true
    }

    if (use_fegenie) {
        ch_fegenie_db = file(params.fegenie_db).exists() ? file(params.fegenie_db) : error("Error: If using --annotate, you must supply prebuilt databases. FEGENIE database file not found at ${params.fegenie_db}")
    }

    if (use_canthyd) {
        ch_canthyd_hmm_db = file(params.canthyd_hmm_db).exists() ? file(params.canthyd_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD HMM database file not found at ${params.canthyd_hmm_db}")
        canthyd_mmseqs_db_path = use_gpu_canthyd ? params.canthyd_mmseqs_gpu_db : params.canthyd_mmseqs_db
        ch_canthyd_mmseqs_db = file(canthyd_mmseqs_db_path).exists() ? file(canthyd_mmseqs_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD MMseqs ${use_gpu_canthyd ? 'gpu' : 'cpu'} database file not found at ${canthyd_mmseqs_db_path}")
        index_mmseqs = true
        ch_canthyd_mmseqs_list = file(params.canthyd_mmseqs_list)
    }

    if (use_vog) {
        ch_vogdb_db = file(params.vog_db).exists() ? file(params.vog_db) : error("Error: If using --annotate, you must supply prebuilt databases. VOG database file not found at ${params.vog_db}")
    }

    if (params.use_viral) {
        viral_db_path = use_gpu_viral ? params.viral_gpu_db : params.viral_db
        ch_viral_db = file(viral_db_path).exists() ? file(viral_db_path) : error("Error: If using --annotate, you must supply prebuilt databases. viral ${use_gpu_viral ? 'gpu' : 'cpu'} database file not found at ${viral_db_path}")
        index_mmseqs = true
    }

    emit:
    ch_kegg_db
    ch_kofam_db
    ch_dbcan_db
    ch_camper_hmm_db
    ch_camper_mmseqs_db
    ch_camper_mmseqs_list
    ch_merops_db
    ch_pfam_mmseqs_db
    ch_heme_db
    ch_sulfur_db
    ch_uniref_db
    ch_metals_db
    ch_antismash_db
    ch_card_db
    ch_card_mmseqs_db
    ch_tcdb_db
    ch_dram_db
    ch_methyl_db
    ch_fegenie_db
    ch_canthyd_hmm_db
    ch_canthyd_mmseqs_db
    ch_canthyd_mmseqs_list
    ch_vogdb_db
    ch_viral_db
    index_mmseqs
}
