//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
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

include { ADD_SQL_DESCRIPTIONS as SQL_UNIREF            } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_VIRAL             } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_MEROPS            } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_KEGG              } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_PFAM              } from "../../modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_DBCAN             } from "../../modules/local/annotate/add_sql_descriptions.nf"

include { HMM_SEARCH as HMM_SEARCH_KOFAM                } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_DBCAN                } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_VOG                  } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CAMPER               } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CANTHYD              } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_SULFUR               } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_FEGENIE              } from "../../modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_METALS               } from "../../modules/local/annotate/hmmsearch.nf"

include { CONCAT_HMM_HITS as CONCAT_HMM_HITS_KOFAM      } from "../../modules/local/annotate/concat_hmm_hits.nf"
include { CONCAT_HMM_HITS as CONCAT_HMM_HITS_VOG        } from "../../modules/local/annotate/concat_hmm_hits.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO DB_SEARCH
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DB_SEARCH {
    take:
    ch_gene_locs  // channel: path(gene_locs_tsv) ]
    ch_called_proteins  // channel: [ val(input_fasta name), path(called_proteins_file.faa) ]
    default_sheet // Path to dummy sheet
    n_fastas // Number of FASTA files to process
    call     // boolean: whether gene calling flag is set
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
    use_vog

    main:

    DB_channel_SETUP(
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
        use_vog

    )

    ch_sql_descriptions_db = file(params.sql_descriptions_db)
    ch_kofam_list = file(params.kofam_list)
    ch_dbcan_fam = file(params.dbcan_fam_activities)
    ch_dbcan_subfam = file(params.dbcan_subfam_activities)
    ch_vog_list = file(params.vog_list)
    ch_camper_hmm_list = file(params.camper_hmm_list)
    ch_canthyd_hmm_list = file(params.cant_hyd_hmm_list)


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


    if (!call) {
        ch_called_proteins = channel
            .fromPath(file(params.input_genes) / params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta file of called genes using --input_genes. Cannot find any called gene fasta files matching: ${params.input_genes}\nNB: Path needs to follow pattern: path/to/directory/" }
            .map { it -> [ it.getBaseName(), it ] }
            
        GENE_LOCS( ch_called_proteins)
        ch_gene_locs = GENE_LOCS.out.prodigal_locs_tsv
        n_fastas = file("$params.input_genes/${params.genes_fmt}").size()
    }

    def formattedOutputchannels = channel.of()

    // Here we will create mmseqs2 index files for each of the inputs if we are going to do a mmseqs2 database
    if (DB_channel_SETUP.out.index_mmseqs) {
        // Use MMSEQS2 to index each called genes protein file
        MMSEQS_INDEX( ch_called_proteins )
        ch_mmseqs_query = MMSEQS_INDEX.out.mmseqs_index_out
    }

    // KEGG annotation
    if (use_kegg) {
        ch_combined_query_locs_kegg = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_KEGG( ch_combined_query_locs_kegg, DB_channel_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, kegg_name )
        ch_kegg_unformatted = MMSEQS_SEARCH_KEGG.out.mmseqs_search_formatted_out

        SQL_KEGG(ch_kegg_unformatted, kegg_name, ch_sql_descriptions_db)
        ch_kegg_formatted = SQL_KEGG.out.sql_formatted_hits

        formattedOutputchannels = formattedOutputchannels.mix(ch_kegg_formatted)
    }
    // KOFAM annotation
    if (use_kofam) {
        if (params.kofam_chunk_size && params.kofam_chunk_size > 0) {
            ch_kofam_chunks = ch_called_proteins
                .splitFasta(by: params.kofam_chunk_size, file: true, elem: 1)
                .combine(ch_gene_locs, by: 0)
            ch_kofam_input = ch_kofam_chunks
                .map { sample, chunk_fasta, gene_locs ->
                    tuple("${sample}___${chunk_fasta.baseName}", chunk_fasta, gene_locs)
                }
            ch_kofam_chunk_to_sample = ch_kofam_chunks
                .map { sample, chunk_fasta, gene_locs ->
                    tuple("${sample}___${chunk_fasta.baseName}", sample)
                }
            HMM_SEARCH_KOFAM (
                ch_kofam_input,
                params.kofam_e_value,
                DB_channel_SETUP.out.ch_kofam_db,
                ch_kofam_list,
                true,
                kofam_name
            )
            ch_kofam_chunk_csvs_per_sample = HMM_SEARCH_KOFAM.out.formatted_hits
                .join(ch_kofam_chunk_to_sample)
                .map { chunk_id, csv, sample -> tuple(sample, csv) }
                .groupTuple()
            CONCAT_HMM_HITS_KOFAM(ch_kofam_chunk_csvs_per_sample, kofam_name)
            ch_kofam_formatted = CONCAT_HMM_HITS_KOFAM.out.combined_hits
        } else {
            ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
            HMM_SEARCH_KOFAM (
                ch_combined_proteins_locs,
                params.kofam_e_value,
                DB_channel_SETUP.out.ch_kofam_db,
                ch_kofam_list,
                true,
                kofam_name
            )
            ch_kofam_formatted = HMM_SEARCH_KOFAM.out.formatted_hits
        }
        formattedOutputchannels = formattedOutputchannels.mix(ch_kofam_formatted)
    }
    // PFAM annotation
    if (use_pfam) {
        ch_combined_query_locs_pfam = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_PFAM( ch_combined_query_locs_pfam, DB_channel_SETUP.out.ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, pfam_name )
        ch_pfam_unformatted = MMSEQS_SEARCH_PFAM.out.mmseqs_search_formatted_out

        SQL_PFAM(ch_pfam_unformatted, pfam_name, ch_sql_descriptions_db)
        ch_pfam_formatted = SQL_PFAM.out.sql_formatted_hits

        formattedOutputchannels = formattedOutputchannels.mix(ch_pfam_formatted)
    }
    // dbCAN annotation
    if  (use_dbcan) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_DBCAN ( 
            ch_combined_proteins_locs, 
            params.dbcan_e_value,
            DB_channel_SETUP.out.ch_dbcan_db,
            default_sheet,
            false,
            dbcan_name
            )
        ch_dbcan_unformatted = HMM_SEARCH_DBCAN.out.formatted_hits
        SQL_DBCAN(ch_dbcan_unformatted, dbcan_name, ch_sql_descriptions_db)
        ch_dbcan_formatted = SQL_DBCAN.out.sql_formatted_hits
        formattedOutputchannels = formattedOutputchannels.mix(ch_dbcan_formatted)
    }
    // CAMPER annotation
    if (use_camper) {
        // HMM
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_CAMPER ( 
            ch_combined_proteins_locs, 
            params.camper_e_value, 
            DB_channel_SETUP.out.ch_camper_hmm_db,
            ch_camper_hmm_list,
            false,
            camper_name
        )
        ch_camper_hmm_formatted = HMM_SEARCH_CAMPER.out.formatted_hits
        formattedOutputchannels = formattedOutputchannels.mix(ch_camper_hmm_formatted)

        // MMseqs
        ch_combined_query_locs_camper = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_CAMPER( ch_combined_query_locs_camper, DB_channel_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_channel_SETUP.out.ch_camper_mmseqs_list, camper_name )
        ch_camper_mmseqs_formatted = MMSEQS_SEARCH_CAMPER.out.mmseqs_search_formatted_out

        formattedOutputchannels = formattedOutputchannels.mix(ch_camper_mmseqs_formatted)
    }
    // FeGenie annotation
    if (use_fegenie) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_FEGENIE ( 
            ch_combined_proteins_locs,  
            params.fegenie_e_value, 
            DB_channel_SETUP.out.ch_fegenie_db,
            default_sheet,
            false,
            fegenie_name
            )
        ch_fegenie_formatted = HMM_SEARCH_FEGENIE.out.formatted_hits
        formattedOutputchannels = formattedOutputchannels.mix(ch_fegenie_formatted)
    }
    // Methyl annotation
    if (use_methyl) {
        ch_combined_query_locs_methyl = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_METHYL( ch_combined_query_locs_methyl, DB_channel_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, methyl_name )
        ch_methyl_mmseqs_formatted = MMSEQS_SEARCH_METHYL.out.mmseqs_search_formatted_out

        formattedOutputchannels = formattedOutputchannels.mix(ch_methyl_mmseqs_formatted)
    }
    // CANT-HYD annotation
    if (use_canthyd) {
        // MMseqs
        ch_combined_query_locs_canthyd = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_CANTHYD( ch_combined_query_locs_canthyd, DB_channel_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_channel_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name )
        ch_canthyd_mmseqs_formatted = MMSEQS_SEARCH_CANTHYD.out.mmseqs_search_formatted_out

        formattedOutputchannels = formattedOutputchannels.mix(ch_canthyd_mmseqs_formatted)

        //HMM
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_CANTHYD ( 
            ch_combined_proteins_locs, 
            params.canthyd_e_value, 
            DB_channel_SETUP.out.ch_canthyd_hmm_db,
            ch_canthyd_hmm_list,
            false,
            canthyd_name
            )
        ch_canthyd_hmm_formatted = HMM_SEARCH_CANTHYD.out.formatted_hits
        formattedOutputchannels = formattedOutputchannels.mix(ch_canthyd_hmm_formatted)

    }
    // Sulfur annotation
    if (use_sulfur) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_SULFUR ( 
            ch_combined_proteins_locs,  
            params.sulfur_e_value, 
            DB_channel_SETUP.out.ch_sulfur_db,
            default_sheet,
            false,
            sulfur_name
            )
        ch_sulfur_formatted = HMM_SEARCH_SULFUR.out.formatted_hits
        formattedOutputchannels = formattedOutputchannels.mix(ch_sulfur_formatted)
    }
    // MEROPS annotation
    if (use_merops) {
        ch_combined_query_locs_merops = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_MEROPS( ch_combined_query_locs_merops, DB_channel_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, merops_name )
        ch_merops_unformatted = MMSEQS_SEARCH_MEROPS.out.mmseqs_search_formatted_out

        SQL_MEROPS(ch_merops_unformatted, merops_name, ch_sql_descriptions_db)
        ch_merops_formatted = SQL_MEROPS.out.sql_formatted_hits

        formattedOutputchannels = formattedOutputchannels.mix(ch_merops_formatted)
    }
    // Uniref annotation
    if (use_uniref) {
        ch_combined_query_locs_uniref = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_UNIREF( ch_combined_query_locs_uniref, DB_channel_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, default_sheet, uniref_name )
        ch_uniref_unformatted = MMSEQS_SEARCH_UNIREF.out.mmseqs_search_formatted_out

        SQL_UNIREF(ch_uniref_unformatted, uniref_name, ch_sql_descriptions_db)
        ch_uniref_formatted = SQL_UNIREF.out.sql_formatted_hits

        formattedOutputchannels = formattedOutputchannels.mix(ch_uniref_formatted)
    }
    // Metals annotation
    if (use_metals) {
        ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
        HMM_SEARCH_METALS ( 
            ch_combined_proteins_locs,  
            params.metals_e_value, 
            DB_channel_SETUP.out.ch_metals_db,
            default_sheet,
            false,
            metals_name
            )
        ch_metals_formatted = HMM_SEARCH_METALS.out.formatted_hits
        formattedOutputchannels = formattedOutputchannels.mix(ch_metals_formatted)
    }
    // VOGdb annotation
    if (use_vog) {
        if (params.vog_chunk_size && params.vog_chunk_size > 0) {
            ch_vog_chunks = ch_called_proteins
                .splitFasta(by: params.vog_chunk_size, file: true, elem: 1)
                .combine(ch_gene_locs, by: 0)
            ch_vog_input = ch_vog_chunks
                .map { sample, chunk_fasta, gene_locs ->
                    tuple("${sample}___${chunk_fasta.baseName}", chunk_fasta, gene_locs)
                }
            ch_vog_chunk_to_sample = ch_vog_chunks
                .map { sample, chunk_fasta, gene_locs ->
                    tuple("${sample}___${chunk_fasta.baseName}", sample)
                }
            HMM_SEARCH_VOG (
                ch_vog_input,
                params.vog_e_value,
                DB_channel_SETUP.out.ch_vogdb_db,
                default_sheet,
                false,
                vogdb_name
            )
            ch_vog_chunk_csvs_per_sample = HMM_SEARCH_VOG.out.formatted_hits
                .join(ch_vog_chunk_to_sample)
                .map { chunk_id, csv, sample -> tuple(sample, csv) }
                .groupTuple()
            CONCAT_HMM_HITS_VOG(ch_vog_chunk_csvs_per_sample, vogdb_name)
            ch_vog_formatted = CONCAT_HMM_HITS_VOG.out.combined_hits
        } else {
            ch_combined_proteins_locs = ch_called_proteins.join(ch_gene_locs)
            HMM_SEARCH_VOG (
                ch_combined_proteins_locs,
                params.vog_e_value,
                DB_channel_SETUP.out.ch_vogdb_db,
                default_sheet,
                false,
                vogdb_name
            )
            ch_vog_formatted = HMM_SEARCH_VOG.out.formatted_hits
        }
        formattedOutputchannels = formattedOutputchannels.mix(ch_vog_formatted)
    }
    // Viral annotation
    if (params.use_viral) {
        ch_combined_query_locs_viral = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_VIRAL( ch_combined_query_locs_viral, DB_channel_SETUP.out.ch_viral_db, params.bit_score_threshold,  params.rbh_bit_score_threshold,default_sheet, viral_name )
        ch_viral_unformatted = MMSEQS_SEARCH_VIRAL.out.mmseqs_search_formatted_out

        SQL_VIRAL(ch_viral_unformatted, viral_name, ch_sql_descriptions_db)
        ch_viral_formatted = SQL_VIRAL.out.sql_formatted_hits

        formattedOutputchannels = formattedOutputchannels.mix(ch_viral_formatted)
    }
    fastas = formattedOutputchannels.map { it -> it[1] }.collect()
    genes = ch_called_proteins.map { it -> it[1] }.collect()
    
    COMBINE_ANNOTATIONS( fastas, genes )
    ch_combined_annotations = COMBINE_ANNOTATIONS.out.combined_annotations_out


    emit:
    ch_combined_annotations  // channel: [ path(combined_annotations_out) ]

}

workflow DB_channel_SETUP {
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
    use_vog


    main:

    index_mmseqs = false
    ch_kegg_db = Channel.empty()
    ch_kofam_db = Channel.empty()
    ch_dbcan_db = Channel.empty()
    ch_camper_hmm_db = Channel.empty()
    ch_camper_mmseqs_db = Channel.empty()
    ch_camper_mmseqs_list = Channel.empty()
    ch_merops_db = Channel.empty()
    ch_pfam_mmseqs_db = Channel.empty()
    ch_heme_db = Channel.empty()
    ch_sulfur_db = Channel.empty()
    ch_uniref_db = Channel.empty()
    ch_metals_db = Channel.empty()
    ch_methyl_db = Channel.empty()
    ch_fegenie_db = Channel.empty()
    ch_canthyd_hmm_db = Channel.empty()
    ch_canthyd_mmseqs_db = Channel.empty()
    ch_canthyd_mmseqs_list = Channel.empty()
    ch_vogdb_db = Channel.empty()
    ch_viral_db = Channel.empty()

    if (use_kegg) {
        ch_kegg_db = file(params.kegg_db).exists() ? file(params.kegg_db) : error("Error: If using --annotate, you must supply prebuilt databases. KEGG database file not found at ${params.kegg_db}")
        index_mmseqs = true
    }

    if (use_kofam) {
        ch_kofam_db = file(params.kofam_db).exists() ? file(params.kofam_db) : error("Error: If using --annotate, you must supply prebuilt databases. KOFAM database file not found at ${params.kofam_db}")
    }

    if (use_dbcan) {
        ch_dbcan_db = file(params.dbcan_db).exists() ? file(params.dbcan_db) : error("Error: If using --annotate, you must supply prebuilt databases. DBCAN database file not found at ${params.dbcan_db}")
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
    ch_methyl_db
    ch_fegenie_db
    ch_canthyd_hmm_db
    ch_canthyd_mmseqs_db
    ch_canthyd_mmseqs_list
    ch_vogdb_db
    ch_viral_db
    index_mmseqs
} 