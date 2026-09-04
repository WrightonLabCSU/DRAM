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
include { MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS_SMALL;  MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS_MEDIUM;  MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS_LARGE  } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL_SMALL;   MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL_MEDIUM;   MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL_LARGE   } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER_SMALL;  MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER_MEDIUM;  MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER_LARGE  } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL_SMALL;  MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL_MEDIUM;  MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL_LARGE  } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD_SMALL; MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD_MEDIUM; MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD_LARGE } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG_SMALL;    MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG_MEDIUM;    MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG_LARGE    } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF_SMALL;  MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF_MEDIUM;  MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF_LARGE  } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM_SMALL;    MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM_MEDIUM;    MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM_LARGE    } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CARD_SMALL;    MMSEQS_SEARCH as MMSEQS_SEARCH_CARD_MEDIUM;    MMSEQS_SEARCH as MMSEQS_SEARCH_CARD_LARGE    } from "../../modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_TCDB_SMALL;    MMSEQS_SEARCH as MMSEQS_SEARCH_TCDB_MEDIUM;    MMSEQS_SEARCH as MMSEQS_SEARCH_TCDB_LARGE    } from "../../modules/local/annotate/mmseqs_search.nf"

include { ADD_SQL_DESCRIPTIONS as SQL_UNIREF            } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_VIRAL             } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_MEROPS            } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_KEGG              } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_PFAM              } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_DBCAN             } from "../../modules/local/annotate/add_sql_descriptions.nf"

include { HMM_SEARCH as HMM_SEARCH_KOFAM_SMALL;   HMM_SEARCH as HMM_SEARCH_KOFAM_MEDIUM;   HMM_SEARCH as HMM_SEARCH_KOFAM_LARGE   } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_DRAM_DB_SMALL; HMM_SEARCH as HMM_SEARCH_DRAM_DB_MEDIUM; HMM_SEARCH as HMM_SEARCH_DRAM_DB_LARGE } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_VOG_SMALL;     HMM_SEARCH as HMM_SEARCH_VOG_MEDIUM;     HMM_SEARCH as HMM_SEARCH_VOG_LARGE     } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CAMPER_SMALL;  HMM_SEARCH as HMM_SEARCH_CAMPER_MEDIUM;  HMM_SEARCH as HMM_SEARCH_CAMPER_LARGE  } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CANTHYD_SMALL; HMM_SEARCH as HMM_SEARCH_CANTHYD_MEDIUM; HMM_SEARCH as HMM_SEARCH_CANTHYD_LARGE } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_SULFUR_SMALL;  HMM_SEARCH as HMM_SEARCH_SULFUR_MEDIUM;  HMM_SEARCH as HMM_SEARCH_SULFUR_LARGE  } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_FEGENIE_SMALL; HMM_SEARCH as HMM_SEARCH_FEGENIE_MEDIUM; HMM_SEARCH as HMM_SEARCH_FEGENIE_LARGE } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_METALS_SMALL;  HMM_SEARCH as HMM_SEARCH_METALS_MEDIUM;  HMM_SEARCH as HMM_SEARCH_METALS_LARGE  } from "../../modules/local/annotate/hmmsearch.nf"

include { ANTISMASH_ANTISMASH                           } from '../../modules/nf-core/antismash/antismash/main'
include { RGI_MAIN                                      } from '../../modules/nf-core/rgi/main/main'
include { RUNDBCAN_EASYSUBSTRATE                        } from '../../modules/nf-core/rundbcan/easysubstrate/main'

include {checkDBVersion                                 } from '../../subworkflows/local/utils_pipeline_setup.nf'
include {resourceBytes; resourceClass                   } from './utils_resource_classes.nf'

def bucketSearchInputs(ch_inputs, database) {
    final long database_bytes = resourceBytes(database)
    ch_inputs
        .map { name, query, gene_locs ->
            tuple(resourceClass(resourceBytes(query) + database_bytes), name, query, gene_locs)
        }
        .branch { resource_class, name, query, gene_locs ->
            small: resource_class == 'small'
            medium: resource_class == 'medium'
            large: resource_class == 'large'
        }
}

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
        ch_kegg_resource = bucketSearchInputs(ch_combined_query_locs_kegg, DB_CHANNEL_SETUP.out.ch_kegg_db.val)
        MMSEQS_SEARCH_KEGG_SMALL(ch_kegg_resource.small, DB_CHANNEL_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, kegg_name)
        MMSEQS_SEARCH_KEGG_MEDIUM(ch_kegg_resource.medium, DB_CHANNEL_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, kegg_name)
        MMSEQS_SEARCH_KEGG_LARGE(ch_kegg_resource.large, DB_CHANNEL_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, kegg_name)
        ch_mmseqs_unformatted = MMSEQS_SEARCH_KEGG_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_KEGG_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_KEGG_LARGE.out.mmseqs_search_formatted_out)

        SQL_KEGG(ch_mmseqs_unformatted, kegg_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_KEGG.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // KOFAM annotation
    if (use_kofam) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_kofam_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_kofam_db.val)
        HMM_SEARCH_KOFAM_SMALL(ch_kofam_resource.small, params.kofam_e_value, DB_CHANNEL_SETUP.out.ch_kofam_db, ch_kofam_list, true, kofam_name)
        HMM_SEARCH_KOFAM_MEDIUM(ch_kofam_resource.medium, params.kofam_e_value, DB_CHANNEL_SETUP.out.ch_kofam_db, ch_kofam_list, true, kofam_name)
        HMM_SEARCH_KOFAM_LARGE(ch_kofam_resource.large, params.kofam_e_value, DB_CHANNEL_SETUP.out.ch_kofam_db, ch_kofam_list, true, kofam_name)
        ch_hmm_formatted = HMM_SEARCH_KOFAM_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_KOFAM_MEDIUM.out.formatted_hits, HMM_SEARCH_KOFAM_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // PFAM annotation
    if (use_pfam) {
        ch_combined_query_locs_pfam = ch_mmseqs_query.join(ch_gene_locs)
        ch_pfam_resource = bucketSearchInputs(ch_combined_query_locs_pfam, DB_CHANNEL_SETUP.out.ch_pfam_mmseqs_db.val)
        MMSEQS_SEARCH_PFAM_SMALL(ch_pfam_resource.small, DB_CHANNEL_SETUP.out.ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, pfam_name)
        MMSEQS_SEARCH_PFAM_MEDIUM(ch_pfam_resource.medium, DB_CHANNEL_SETUP.out.ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, pfam_name)
        MMSEQS_SEARCH_PFAM_LARGE(ch_pfam_resource.large, DB_CHANNEL_SETUP.out.ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, pfam_name)
        ch_mmseqs_unformatted = MMSEQS_SEARCH_PFAM_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_PFAM_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_PFAM_LARGE.out.mmseqs_search_formatted_out)

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
        ch_camper_hmm_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_camper_hmm_db.val)
        HMM_SEARCH_CAMPER_SMALL(ch_camper_hmm_resource.small, params.camper_e_value, DB_CHANNEL_SETUP.out.ch_camper_hmm_db, ch_camper_hmm_list, false, camper_name)
        HMM_SEARCH_CAMPER_MEDIUM(ch_camper_hmm_resource.medium, params.camper_e_value, DB_CHANNEL_SETUP.out.ch_camper_hmm_db, ch_camper_hmm_list, false, camper_name)
        HMM_SEARCH_CAMPER_LARGE(ch_camper_hmm_resource.large, params.camper_e_value, DB_CHANNEL_SETUP.out.ch_camper_hmm_db, ch_camper_hmm_list, false, camper_name)
        ch_hmm_formatted = HMM_SEARCH_CAMPER_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_CAMPER_MEDIUM.out.formatted_hits, HMM_SEARCH_CAMPER_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)

        // MMseqs
        ch_combined_query_locs_camper = ch_mmseqs_query.join(ch_gene_locs)
        ch_camper_mmseqs_resource = bucketSearchInputs(ch_combined_query_locs_camper, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db.val)
        MMSEQS_SEARCH_CAMPER_SMALL(ch_camper_mmseqs_resource.small, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_list, camper_name)
        MMSEQS_SEARCH_CAMPER_MEDIUM(ch_camper_mmseqs_resource.medium, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_list, camper_name)
        MMSEQS_SEARCH_CAMPER_LARGE(ch_camper_mmseqs_resource.large, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_list, camper_name)
        ch_mmseqs_formatted = MMSEQS_SEARCH_CAMPER_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_CAMPER_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_CAMPER_LARGE.out.mmseqs_search_formatted_out)
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // FeGenie annotation
    if (use_fegenie) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_fegenie_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_fegenie_db.val)
        HMM_SEARCH_FEGENIE_SMALL(ch_fegenie_resource.small, params.fegenie_e_value, DB_CHANNEL_SETUP.out.ch_fegenie_db, default_sheet, false, fegenie_name)
        HMM_SEARCH_FEGENIE_MEDIUM(ch_fegenie_resource.medium, params.fegenie_e_value, DB_CHANNEL_SETUP.out.ch_fegenie_db, default_sheet, false, fegenie_name)
        HMM_SEARCH_FEGENIE_LARGE(ch_fegenie_resource.large, params.fegenie_e_value, DB_CHANNEL_SETUP.out.ch_fegenie_db, default_sheet, false, fegenie_name)
        ch_hmm_formatted = HMM_SEARCH_FEGENIE_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_FEGENIE_MEDIUM.out.formatted_hits, HMM_SEARCH_FEGENIE_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // Methyl annotation
    if (use_methyl) {
        ch_combined_query_locs_methyl = ch_mmseqs_query.join(ch_gene_locs)
        ch_methyl_resource = bucketSearchInputs(ch_combined_query_locs_methyl, DB_CHANNEL_SETUP.out.ch_methyl_db.val)
        MMSEQS_SEARCH_METHYL_SMALL(ch_methyl_resource.small, DB_CHANNEL_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, methyl_name)
        MMSEQS_SEARCH_METHYL_MEDIUM(ch_methyl_resource.medium, DB_CHANNEL_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, methyl_name)
        MMSEQS_SEARCH_METHYL_LARGE(ch_methyl_resource.large, DB_CHANNEL_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, methyl_name)
        ch_mmseqs_formatted = MMSEQS_SEARCH_METHYL_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_METHYL_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_METHYL_LARGE.out.mmseqs_search_formatted_out)
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // CANT-HYD annotation
    if (use_canthyd) {
        // MMseqs
        ch_combined_query_locs_canthyd = ch_mmseqs_query.join(ch_gene_locs)
        ch_canthyd_mmseqs_resource = bucketSearchInputs(ch_combined_query_locs_canthyd, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db.val)
        MMSEQS_SEARCH_CANTHYD_SMALL(ch_canthyd_mmseqs_resource.small, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name)
        MMSEQS_SEARCH_CANTHYD_MEDIUM(ch_canthyd_mmseqs_resource.medium, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name)
        MMSEQS_SEARCH_CANTHYD_LARGE(ch_canthyd_mmseqs_resource.large, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name)
        ch_mmseqs_formatted = MMSEQS_SEARCH_CANTHYD_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_CANTHYD_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_CANTHYD_LARGE.out.mmseqs_search_formatted_out)
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)

        //HMM
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_canthyd_hmm_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_canthyd_hmm_db.val)
        HMM_SEARCH_CANTHYD_SMALL(ch_canthyd_hmm_resource.small, params.canthyd_e_value, DB_CHANNEL_SETUP.out.ch_canthyd_hmm_db, ch_canthyd_hmm_list, false, canthyd_name)
        HMM_SEARCH_CANTHYD_MEDIUM(ch_canthyd_hmm_resource.medium, params.canthyd_e_value, DB_CHANNEL_SETUP.out.ch_canthyd_hmm_db, ch_canthyd_hmm_list, false, canthyd_name)
        HMM_SEARCH_CANTHYD_LARGE(ch_canthyd_hmm_resource.large, params.canthyd_e_value, DB_CHANNEL_SETUP.out.ch_canthyd_hmm_db, ch_canthyd_hmm_list, false, canthyd_name)
        ch_hmm_formatted = HMM_SEARCH_CANTHYD_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_CANTHYD_MEDIUM.out.formatted_hits, HMM_SEARCH_CANTHYD_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // Sulfur annotation
    if (use_sulfur) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_sulfur_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_sulfur_db.val)
        HMM_SEARCH_SULFUR_SMALL(ch_sulfur_resource.small, params.sulfur_e_value, DB_CHANNEL_SETUP.out.ch_sulfur_db, default_sheet, false, sulfur_name)
        HMM_SEARCH_SULFUR_MEDIUM(ch_sulfur_resource.medium, params.sulfur_e_value, DB_CHANNEL_SETUP.out.ch_sulfur_db, default_sheet, false, sulfur_name)
        HMM_SEARCH_SULFUR_LARGE(ch_sulfur_resource.large, params.sulfur_e_value, DB_CHANNEL_SETUP.out.ch_sulfur_db, default_sheet, false, sulfur_name)
        ch_hmm_formatted = HMM_SEARCH_SULFUR_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_SULFUR_MEDIUM.out.formatted_hits, HMM_SEARCH_SULFUR_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // MEROPS annotation
    if (use_merops) {
        ch_combined_query_locs_merops = ch_mmseqs_query.join(ch_gene_locs)
        ch_merops_resource = bucketSearchInputs(ch_combined_query_locs_merops, DB_CHANNEL_SETUP.out.ch_merops_db.val)
        MMSEQS_SEARCH_MEROPS_SMALL(ch_merops_resource.small, DB_CHANNEL_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, merops_name)
        MMSEQS_SEARCH_MEROPS_MEDIUM(ch_merops_resource.medium, DB_CHANNEL_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, merops_name)
        MMSEQS_SEARCH_MEROPS_LARGE(ch_merops_resource.large, DB_CHANNEL_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, merops_name)
        ch_mmseqs_unformatted = MMSEQS_SEARCH_MEROPS_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_MEROPS_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_MEROPS_LARGE.out.mmseqs_search_formatted_out)

        SQL_MEROPS(ch_mmseqs_unformatted, merops_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_MEROPS.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // Uniref annotation
    if (use_uniref) {
        ch_combined_query_locs_uniref = ch_mmseqs_query.join(ch_gene_locs)
        ch_uniref_resource = bucketSearchInputs(ch_combined_query_locs_uniref, DB_CHANNEL_SETUP.out.ch_uniref_db.val)
        MMSEQS_SEARCH_UNIREF_SMALL(ch_uniref_resource.small, DB_CHANNEL_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, uniref_name)
        MMSEQS_SEARCH_UNIREF_MEDIUM(ch_uniref_resource.medium, DB_CHANNEL_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, uniref_name)
        MMSEQS_SEARCH_UNIREF_LARGE(ch_uniref_resource.large, DB_CHANNEL_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, uniref_name)
        ch_mmseqs_unformatted = MMSEQS_SEARCH_UNIREF_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_UNIREF_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_UNIREF_LARGE.out.mmseqs_search_formatted_out)

        SQL_UNIREF(ch_mmseqs_unformatted, uniref_name, ch_sql_descriptions_db)
        ch_mmseqs_formatted = SQL_UNIREF.out.sql_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // Metals annotation
    if (use_metals) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_metals_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_metals_db.val)
        HMM_SEARCH_METALS_SMALL(ch_metals_resource.small, params.metals_e_value, DB_CHANNEL_SETUP.out.ch_metals_db, default_sheet, false, metals_name)
        HMM_SEARCH_METALS_MEDIUM(ch_metals_resource.medium, params.metals_e_value, DB_CHANNEL_SETUP.out.ch_metals_db, default_sheet, false, metals_name)
        HMM_SEARCH_METALS_LARGE(ch_metals_resource.large, params.metals_e_value, DB_CHANNEL_SETUP.out.ch_metals_db, default_sheet, false, metals_name)
        ch_hmm_formatted = HMM_SEARCH_METALS_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_METALS_MEDIUM.out.formatted_hits, HMM_SEARCH_METALS_LARGE.out.formatted_hits)
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
        ch_card_resource = bucketSearchInputs(ch_combined_query_locs_card, DB_CHANNEL_SETUP.out.ch_card_db.val)
        MMSEQS_SEARCH_CARD_SMALL(ch_card_resource.small, DB_CHANNEL_SETUP.out.ch_card_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, card_name)
        MMSEQS_SEARCH_CARD_MEDIUM(ch_card_resource.medium, DB_CHANNEL_SETUP.out.ch_card_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, card_name)
        MMSEQS_SEARCH_CARD_LARGE(ch_card_resource.large, DB_CHANNEL_SETUP.out.ch_card_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, card_name)
        ch_mmseqs_formatted = MMSEQS_SEARCH_CARD_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_CARD_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_CARD_LARGE.out.mmseqs_search_formatted_out)
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // TCDB annotation
    if (use_tcdb) {
        ch_combined_query_locs_tcdb = ch_mmseqs_query.join(ch_gene_locs)
        ch_tcdb_resource = bucketSearchInputs(ch_combined_query_locs_tcdb, DB_CHANNEL_SETUP.out.ch_tcdb_db.val)
        MMSEQS_SEARCH_TCDB_SMALL(ch_tcdb_resource.small, DB_CHANNEL_SETUP.out.ch_tcdb_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, tcdb_name)
        MMSEQS_SEARCH_TCDB_MEDIUM(ch_tcdb_resource.medium, DB_CHANNEL_SETUP.out.ch_tcdb_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, tcdb_name)
        MMSEQS_SEARCH_TCDB_LARGE(ch_tcdb_resource.large, DB_CHANNEL_SETUP.out.ch_tcdb_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, tcdb_name)
        ch_mmseqs_formatted = MMSEQS_SEARCH_TCDB_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_TCDB_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_TCDB_LARGE.out.mmseqs_search_formatted_out)
        formattedOutputChannels = formattedOutputChannels.mix(ch_mmseqs_formatted)
    }
    // DRAM DB annotation
    if (use_dram_db) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_dram_db_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_dram_db.val)
        HMM_SEARCH_DRAM_DB_SMALL(ch_dram_db_resource.small, "", DB_CHANNEL_SETUP.out.ch_dram_db, ch_dram_db_hmm_list, false, dram_db_name)
        HMM_SEARCH_DRAM_DB_MEDIUM(ch_dram_db_resource.medium, "", DB_CHANNEL_SETUP.out.ch_dram_db, ch_dram_db_hmm_list, false, dram_db_name)
        HMM_SEARCH_DRAM_DB_LARGE(ch_dram_db_resource.large, "", DB_CHANNEL_SETUP.out.ch_dram_db, ch_dram_db_hmm_list, false, dram_db_name)
        ch_hmm_formatted = HMM_SEARCH_DRAM_DB_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_DRAM_DB_MEDIUM.out.formatted_hits, HMM_SEARCH_DRAM_DB_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // VOGdb annotation
    if (use_vog) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        ch_vog_resource = bucketSearchInputs(ch_combined_proteins_locs, DB_CHANNEL_SETUP.out.ch_vogdb_db.val)
        HMM_SEARCH_VOG_SMALL(ch_vog_resource.small, params.vog_e_value, DB_CHANNEL_SETUP.out.ch_vogdb_db, default_sheet, false, vogdb_name)
        HMM_SEARCH_VOG_MEDIUM(ch_vog_resource.medium, params.vog_e_value, DB_CHANNEL_SETUP.out.ch_vogdb_db, default_sheet, false, vogdb_name)
        HMM_SEARCH_VOG_LARGE(ch_vog_resource.large, params.vog_e_value, DB_CHANNEL_SETUP.out.ch_vogdb_db, default_sheet, false, vogdb_name)
        ch_hmm_formatted = HMM_SEARCH_VOG_SMALL.out.formatted_hits
            .mix(HMM_SEARCH_VOG_MEDIUM.out.formatted_hits, HMM_SEARCH_VOG_LARGE.out.formatted_hits)
        formattedOutputChannels = formattedOutputChannels.mix(ch_hmm_formatted)
    }
    // Viral annotation
    if (params.use_viral) {
        ch_combined_query_locs_viral = ch_mmseqs_query.join(ch_gene_locs)
        ch_viral_resource = bucketSearchInputs(ch_combined_query_locs_viral, DB_CHANNEL_SETUP.out.ch_viral_db.val)
        MMSEQS_SEARCH_VIRAL_SMALL(ch_viral_resource.small, DB_CHANNEL_SETUP.out.ch_viral_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, viral_name)
        MMSEQS_SEARCH_VIRAL_MEDIUM(ch_viral_resource.medium, DB_CHANNEL_SETUP.out.ch_viral_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, viral_name)
        MMSEQS_SEARCH_VIRAL_LARGE(ch_viral_resource.large, DB_CHANNEL_SETUP.out.ch_viral_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, viral_name)
        ch_mmseqs_unformatted = MMSEQS_SEARCH_VIRAL_SMALL.out.mmseqs_search_formatted_out
            .mix(MMSEQS_SEARCH_VIRAL_MEDIUM.out.mmseqs_search_formatted_out, MMSEQS_SEARCH_VIRAL_LARGE.out.mmseqs_search_formatted_out)

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
        ch_kegg_db = file(params.kegg_db).exists() ? file(params.kegg_db) : error("Error: If using --annotate, you must supply prebuilt databases. KEGG database file not found at ${params.kegg_db}")
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
        ch_camper_mmseqs_db = file(params.camper_mmseqs_db).exists() ? file(params.camper_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER MMseqs2 database file not found at ${params.camper_mmseqs_db}")
        index_mmseqs = true
        ch_camper_mmseqs_list = file(params.camper_mmseqs_list)
    }

    if (use_merops) {
        ch_merops_db = file(params.merops_db).exists() ? file(params.merops_db) : error("Error: If using --annotate, you must supply prebuilt databases. MEROPS database file not found at ${params.merops_db}")
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
        ch_uniref_db = file(params.uniref_db).exists() ? file(params.uniref_db) : error("Error: If using --annotate, you must supply prebuilt databases. UNIREF database file not found at ${params.uniref_db}")
        index_mmseqs = true
    }

    if (use_metals) {
        ch_metals_db = file(params.metals_db).exists() ? file(params.metals_db) : error("Error: If using --annotate, you must supply prebuilt databases. METALS database file not found at ${params.metals_db}")
    }

    if (use_antismash) {
        ch_antismash_db = file(params.antismash_db).exists() ? file(params.antismash_db) : error("Error: If using --annotate, you must supply prebuilt databases. antismash database file not found at ${params.antismash_db}")
    }

    if (use_rgi || use_card) {
        ch_card_db = file(params.card_db).exists() ? file(params.card_db) : error("Error: If using --annotate, you must supply prebuilt databases. rgi database file not found at ${params.card_db}")
        // rgi software uses the raw fasta, but card search we use the mmseqs database
        if (use_card) {
            index_mmseqs = true
        }
    }

    if (use_tcdb) {
        ch_tcdb_db = file(params.tcdb_db).exists() ? file(params.tcdb_db) : error("Error: If using --annotate, you must supply prebuilt databases. tcdb database file not found at ${params.tcdb_db}")
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
        ch_methyl_db = file(params.methyl_db).exists() ? file(params.methyl_db) : error("Error: If using --annotate, you must supply prebuilt databases. METHYL database file not found at ${params.methyl_db}")
        index_mmseqs = true
    }

    if (use_fegenie) {
        ch_fegenie_db = file(params.fegenie_db).exists() ? file(params.fegenie_db) : error("Error: If using --annotate, you must supply prebuilt databases. FEGENIE database file not found at ${params.fegenie_db}")
    }

    if (use_canthyd) {
        ch_canthyd_hmm_db = file(params.canthyd_hmm_db).exists() ? file(params.canthyd_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD HMM database file not found at ${params.canthyd_hmm_db}")
        ch_canthyd_mmseqs_db = file(params.canthyd_mmseqs_db).exists() ? file(params.canthyd_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD MMseqs database file not found at ${params.canthyd_mmseqs_db}")
        index_mmseqs = true
        ch_canthyd_mmseqs_list = file(params.canthyd_mmseqs_list)
    }

    if (use_vog) {
        ch_vogdb_db = file(params.vog_db).exists() ? file(params.vog_db) : error("Error: If using --annotate, you must supply prebuilt databases. VOG database file not found at ${params.vog_db}")
    }

    if (params.use_viral) {
        ch_viral_db = file(params.viral_db).exists() ? file(params.viral_db) : error("Error: If using --annotate, you must supply prebuilt databases. viral database file not found at ${params.viral_db}")
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
