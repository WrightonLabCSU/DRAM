//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENE_LOCS                                     } from "${projectDir}/modules/local/annotate/gene_locs.nf"

include { GENERIC_HMM_FORMATTER                         } from "${projectDir}/modules/local/annotate/generic_hmm_formatter.nf"  // TODO, This has hard coded paths on the server to the python file. Need to fix this.
include { KEGG_HMM_FORMATTER                            } from "${projectDir}/modules/local/annotate/kegg_hmm_formatter.nf"
include { KOFAM_HMM_FORMATTER                           } from "${projectDir}/modules/local/annotate/kofam_hmm_formatter.nf"
include { DBCAN_HMM_FORMATTER                           } from "${projectDir}/modules/local/annotate/dbcan_hmm_formatter.nf"
include { VOG_HMM_FORMATTER                             } from "${projectDir}/modules/local/annotate/vog_hmm_formatter.nf"
include { CAMPER_HMM_FORMATTER                          } from "${projectDir}/modules/local/annotate/camper_hmm_formatter.nf"
include { CANTHYD_HMM_FORMATTER                         } from "${projectDir}/modules/local/annotate/canthyd_hmm_formatter.nf"
include { SULFUR_HMM_FORMATTER                          } from "${projectDir}/modules/local/annotate/sulfur_hmm_formatter.nf"
include { FEGENIE_HMM_FORMATTER                         } from "${projectDir}/modules/local/annotate/fegenie_hmm_formatter.nf"

include { COMBINE_ANNOTATIONS                           } from "${projectDir}/modules/local/annotate/combine_annotations.nf"

include { MMSEQS_INDEX                                  } from "${projectDir}/modules/local/annotate/mmseqs_index.nf"

// NextFlow only process with the same name in the same workflow, so either alias it or include it a different workflow
include { MMSEQS_SEARCH as MMSEQS_SEARCH_MEROPS         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_VIRAL          } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CAMPER         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_METHYL         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_CANTHYD        } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_KEGG           } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_UNIREF         } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"
include { MMSEQS_SEARCH as MMSEQS_SEARCH_PFAM           } from "${projectDir}/modules/local/annotate/mmseqs_search.nf"

include { ADD_SQL_DESCRIPTIONS as SQL_UNIREF            } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_VIRAL             } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_MEROPS            } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_KEGG              } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"
include { ADD_SQL_DESCRIPTIONS as SQL_PFAM              } from "${projectDir}/modules/local/annotate/add_sql_descriptions.nf"

include { HMM_SEARCH as HMM_SEARCH_KOFAM                } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_DBCAN                } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_VOG                  } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CAMPER               } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_CANTHYD              } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_SULFUR               } from "${projectDir}/modules/local/annotate/hmmsearch.nf"
include { HMM_SEARCH as HMM_SEARCH_FEGENIE              } from "${projectDir}/modules/local/annotate/hmmsearch.nf"

include { PARSE_HMM as PARSE_HMM_KOFAM                  } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_DBCAN                  } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_VOG                    } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_CAMPER                 } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_CANTHYD                } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_SULFUR                 } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"
include { PARSE_HMM as PARSE_HMM_FEGENIE                } from "${projectDir}/modules/local/annotate/parse_hmmsearch.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO ANNOTATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    ch_fasta  // channel: [ val(sample name), path(fasta) ]
    ch_gene_locs  // channel: [ val(sample name), path(gene_locs_tsv) ]
    ch_called_proteins  // channel: [ val(sample name), path(called_proteins_file.faa) ]
    ch_dummy_sheet // Path to dummy sheet

    main:

    DB_CHANNEL_SETUP()

    ch_sql_descriptions_db = file(params.sql_descriptions_db)
    ch_kofam_list = file(params.kofam_list)
    ch_canthyd_list = file(params.cant_hyd_hmm_list)
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


    if (!params.call) {
        ch_called_proteins = Channel
            .fromPath(file(params.input_genes) / params.genes_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "If you specify --annotate without --call, you must provide a fasta file of called genes using --input_genes. Cannot find any called gene fasta files matching: ${params.input_genes}\nNB: Path needs to follow pattern: path/to/directory/" }
            .map {
                sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
                tuple(sampleName, it)
            }

        GENE_LOCS( ch_called_proteins)
        ch_gene_locs = GENE_LOCS.out.prodigal_locs_tsv
    }

    def formattedOutputChannels = channel.of()

    // Here we will create mmseqs2 index files for each of the inputs if we are going to do a mmseqs2 database
    if (DB_CHANNEL_SETUP.out.index_mmseqs) {
        // Use MMSEQS2 to index each called genes protein file
        MMSEQS_INDEX( ch_called_proteins )
        ch_mmseqs_query = MMSEQS_INDEX.out.mmseqs_index_out
    }

    // KEGG annotation
    if (params.use_kegg) {
        ch_combined_query_locs_kegg = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_KEGG( ch_combined_query_locs_kegg, DB_CHANNEL_SETUP.out.ch_kegg_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, kegg_name )
        ch_kegg_unformatted = MMSEQS_SEARCH_KEGG.out.mmseqs_search_formatted_out

        SQL_KEGG(ch_kegg_unformatted, kegg_name, ch_sql_descriptions_db)
        ch_kegg_formatted = SQL_KEGG.out.sql_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_kegg_formatted)
    }
    // KOFAM annotation
    if (params.use_kofam) {
        HMM_SEARCH_KOFAM ( ch_called_proteins, params.kofam_e_value, DB_CHANNEL_SETUP.out.ch_kofam_db )
        ch_kofam_hmms = HMM_SEARCH_KOFAM.out.hmm_search_out

        PARSE_HMM_KOFAM ( ch_kofam_hmms )
        ch_kofam_parsed = PARSE_HMM_KOFAM.out.parsed_hmm

        ch_combined_hits_locs_kofam = ch_kofam_parsed.join(ch_gene_locs)
        KOFAM_HMM_FORMATTER ( ch_combined_hits_locs_kofam, ch_kofam_list )
        ch_kofam_formatted = KOFAM_HMM_FORMATTER.out.kofam_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_kofam_formatted)
    }
    // PFAM annotation
    if (params.use_pfam) {
        ch_combined_query_locs_pfam = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_PFAM( ch_combined_query_locs_pfam, DB_CHANNEL_SETUP.out.ch_pfam_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, pfam_name )
        ch_pfam_unformatted = MMSEQS_SEARCH_PFAM.out.mmseqs_search_formatted_out

        SQL_PFAM(ch_pfam_unformatted, pfam_name, ch_sql_descriptions_db)
        ch_pfam_formatted = SQL_PFAM.out.sql_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_pfam_formatted)
    }
    // dbCAN annotation
    if  (params.use_dbcan) {

        HMM_SEARCH_DBCAN ( ch_called_proteins, params.dbcan_e_value , DB_CHANNEL_SETUP.out.ch_dbcan_db)
        ch_dbcan_hmms = HMM_SEARCH_DBCAN.out.hmm_search_out

        PARSE_HMM_DBCAN ( ch_dbcan_hmms )
        ch_dbcan_parsed = PARSE_HMM_DBCAN.out.parsed_hmm

        ch_combined_hits_locs_dbcan = ch_dbcan_parsed.join(ch_gene_locs)
        DBCAN_HMM_FORMATTER ( ch_combined_hits_locs_dbcan, dbcan_name, ch_sql_descriptions_db )
        ch_dbcan_formatted = DBCAN_HMM_FORMATTER.out.sql_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_dbcan_formatted)
    }
    // CAMPER annotation
    if (params.use_camper) {
        // HMM
        HMM_SEARCH_CAMPER ( ch_called_proteins, params.camper_e_value , DB_CHANNEL_SETUP.out.ch_camper_hmm_db)
        ch_camper_hmms = HMM_SEARCH_CAMPER.out.hmm_search_out

        PARSE_HMM_CAMPER ( ch_camper_hmms )
        ch_camper_parsed = PARSE_HMM_CAMPER.out.parsed_hmm

        ch_combined_hits_locs_camper = ch_camper_parsed.join(ch_gene_locs)
        CAMPER_HMM_FORMATTER ( ch_combined_hits_locs_camper, ch_camper_hmm_list )
        ch_camper_hmm_formatted = CAMPER_HMM_FORMATTER.out.camper_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_camper_hmm_formatted)

        // MMseqs
        ch_combined_query_locs_camper = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_CAMPER( ch_combined_query_locs_camper, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_camper_mmseqs_list, camper_name )
        ch_camper_mmseqs_formatted = MMSEQS_SEARCH_CAMPER.out.mmseqs_search_formatted_out

        formattedOutputChannels = formattedOutputChannels.mix(ch_camper_mmseqs_formatted)
    }
    // FeGenie annotation
    if (params.use_fegenie) {
        HMM_SEARCH_FEGENIE ( ch_called_proteins,  params.fegenie_e_value, DB_CHANNEL_SETUP.out.ch_fegenie_db )
        ch_fegenie_hmms = HMM_SEARCH_FEGENIE.out.hmm_search_out

        PARSE_HMM_FEGENIE ( ch_fegenie_hmms )
        ch_fegenie_parsed = PARSE_HMM_FEGENIE.out.parsed_hmm

        ch_combined_hits_locs_fegenie = ch_fegenie_parsed.join(ch_gene_locs)
        FEGENIE_HMM_FORMATTER ( ch_combined_hits_locs_fegenie )
        ch_fegenie_formatted = FEGENIE_HMM_FORMATTER.out.fegenie_formatted_hits
        formattedOutputChannels = formattedOutputChannels.mix(ch_fegenie_formatted)
    }
    // Methyl annotation
    if (params.use_methyl) {
        ch_combined_query_locs_methyl = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_METHYL( ch_combined_query_locs_methyl, DB_CHANNEL_SETUP.out.ch_methyl_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, methyl_name )
        ch_methyl_mmseqs_formatted = MMSEQS_SEARCH_METHYL.out.mmseqs_search_formatted_out

        formattedOutputChannels = formattedOutputChannels.mix(ch_methyl_mmseqs_formatted)
    }
    // CANT-HYD annotation
    if (params.use_canthyd) {
        // MMseqs
        ch_combined_query_locs_canthyd = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_CANTHYD( ch_combined_query_locs_canthyd, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_db, params.bit_score_threshold, params.rbh_bit_score_threshold, DB_CHANNEL_SETUP.out.ch_canthyd_mmseqs_list, canthyd_name )
        ch_canthyd_mmseqs_formatted = MMSEQS_SEARCH_CANTHYD.out.mmseqs_search_formatted_out

        formattedOutputChannels = formattedOutputChannels.mix(ch_canthyd_mmseqs_formatted)

        //HMM
        HMM_SEARCH_CANTHYD ( ch_called_proteins, params.canthyd_e_value , DB_CHANNEL_SETUP.out.ch_canthyd_hmm_db)
        ch_canthyd_hmms = HMM_SEARCH_CANTHYD.out.hmm_search_out

        PARSE_HMM_CANTHYD ( ch_canthyd_hmms )
        ch_canthyd_parsed = PARSE_HMM_CANTHYD.out.parsed_hmm

        ch_combined_hits_locs_canthyd = ch_canthyd_parsed.join(ch_gene_locs)
        CANTHYD_HMM_FORMATTER ( ch_combined_hits_locs_canthyd, ch_canthyd_hmm_list )
        ch_canthyd_hmm_formatted = CANTHYD_HMM_FORMATTER.out.canthyd_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_canthyd_hmm_formatted)

    }
    // Sulfur annotation
    if (params.use_sulfur) {
        HMM_SEARCH_SULFUR ( ch_called_proteins,  params.sulfur_e_value, DB_CHANNEL_SETUP.out.ch_sulfur_db )
        ch_sulfur_hmms = HMM_SEARCH_SULFUR.out.hmm_search_out

        PARSE_HMM_SULFUR ( ch_sulfur_hmms )
        ch_sulfur_parsed = PARSE_HMM_SULFUR.out.parsed_hmm

        ch_combined_hits_locs_sulfur = ch_sulfur_parsed.join(ch_gene_locs)
        SULFUR_HMM_FORMATTER ( ch_combined_hits_locs_sulfur )
        ch_sulfur_formatted = SULFUR_HMM_FORMATTER.out.sulfur_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_sulfur_formatted)
    }
    // MEROPS annotation
    if (params.use_merops) {
        ch_combined_query_locs_merops = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_MEROPS( ch_combined_query_locs_merops, DB_CHANNEL_SETUP.out.ch_merops_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, merops_name )
        ch_merops_unformatted = MMSEQS_SEARCH_MEROPS.out.mmseqs_search_formatted_out

        SQL_MEROPS(ch_merops_unformatted, merops_name, ch_sql_descriptions_db)
        ch_merops_formatted = SQL_MEROPS.out.sql_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_merops_formatted)
    }
    // Uniref annotation
    if (params.use_uniref) {
        ch_combined_query_locs_uniref = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_UNIREF( ch_combined_query_locs_uniref, DB_CHANNEL_SETUP.out.ch_uniref_db, params.bit_score_threshold, params.rbh_bit_score_threshold, ch_dummy_sheet, uniref_name )
        ch_uniref_unformatted = MMSEQS_SEARCH_UNIREF.out.mmseqs_search_formatted_out

        SQL_UNIREF(ch_uniref_unformatted, uniref_name, ch_sql_descriptions_db)
        ch_uniref_formatted = SQL_UNIREF.out.sql_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_uniref_formatted)
    }
    // VOGdb annotation
    if (params.use_vogdb) {
        HMM_SEARCH_VOG ( ch_called_proteins, params.vog_e_value , DB_CHANNEL_SETUP.out.ch_vogdb_db )
        ch_vog_hmms = HMM_SEARCH_VOG.out.hmm_search_out

        PARSE_HMM_VOG ( ch_vog_hmms )
        ch_vog_parsed = PARSE_HMM_VOG.out.parsed_hmm

        ch_combined_hits_locs_vog = ch_vog_parsed.join(ch_gene_locs)
        VOG_HMM_FORMATTER ( ch_combined_hits_locs_vog, vogdb_name, ch_sql_descriptions_db )
        ch_vog_formatted = VOG_HMM_FORMATTER.out.vog_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_vog_formatted)
    }
    // Viral annotation
    if (params.use_viral) {
        ch_combined_query_locs_viral = ch_mmseqs_query.join(ch_gene_locs)
        MMSEQS_SEARCH_VIRAL( ch_combined_query_locs_viral, DB_CHANNEL_SETUP.out.ch_viral_db, params.bit_score_threshold,  params.rbh_bit_score_threshold,ch_dummy_sheet, viral_name )
        ch_viral_unformatted = MMSEQS_SEARCH_VIRAL.out.mmseqs_search_formatted_out

        SQL_VIRAL(ch_viral_unformatted, viral_name, ch_sql_descriptions_db)
        ch_viral_formatted = SQL_VIRAL.out.sql_formatted_hits

        formattedOutputChannels = formattedOutputChannels.mix(ch_viral_formatted)
    }


    // Collect all formatted annotation output files
    Channel.empty()
        .mix( formattedOutputChannels )
        .collect()
        .set { collected_formatted_hits }

    // COMBINE_ANNOTATIONS collects all annotations files across ALL databases
    COMBINE_ANNOTATIONS( collected_formatted_hits )
    ch_combined_annotations = COMBINE_ANNOTATIONS.out.combined_annotations_out


    emit:
    ch_combined_annotations  // channel: [ path(combined_annotations_out) ]

}

workflow DB_CHANNEL_SETUP {
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
    ch_methyl_db = Channel.empty()
    ch_fegenie_db = Channel.empty()
    ch_canthyd_hmm_db = Channel.empty()
    ch_canthyd_mmseqs_db = Channel.empty()
    ch_canthyd_mmseqs_list = Channel.empty()
    ch_vogdb_db = Channel.empty()
    ch_viral_db = Channel.empty()

    if (params.use_kegg) {
        ch_kegg_db = file(params.kegg_db).exists() ? file(params.kegg_db) : error("Error: If using --annotate, you must supply prebuilt databases. KEGG database file not found at ${params.kegg_db}")
        index_mmseqs = true
    }

    if (params.use_kofam) {
        ch_kofam_db = file(params.kofam_db).exists() ? file(params.kofam_db) : error("Error: If using --annotate, you must supply prebuilt databases. KOFAM database file not found at ${params.kofam_db}")
    }

    if (params.use_dbcan) {
        ch_dbcan_db = file(params.dbcan_db).exists() ? file(params.dbcan_db) : error("Error: If using --annotate, you must supply prebuilt databases. DBCAN database file not found at ${params.dbcan_db}")
    }

    if (params.use_camper) {
        ch_camper_hmm_db = file(params.camper_hmm_db).exists() ? file(params.camper_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER HMM database file not found at ${params.camper_hmm_db}")
        ch_camper_mmseqs_db = file(params.camper_mmseqs_db).exists() ? file(params.camper_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CAMPER MMseqs2 database file not found at ${params.camper_mmseqs_db}")
        index_mmseqs = true
        ch_camper_mmseqs_list = file(params.camper_mmseqs_list)
    }

    if (params.use_merops) {
        ch_merops_db = file(params.merops_db).exists() ? file(params.merops_db) : error("Error: If using --annotate, you must supply prebuilt databases. MEROPS database file not found at ${params.merops_db}")
        index_mmseqs = true
    }

    if (params.use_pfam) {
        ch_pfam_mmseqs_db = file(params.pfam_mmseq_db).exists() ? file(params.pfam_mmseq_db) : error("Error: If using --annotate, you must supply prebuilt databases. PFAM database file not found at ${params.pfam_mmseq_db}")
        index_mmseqs = true
    }

    // if (params.use_heme) {
    //     ch_heme_db = file(params.heme_db).exists() ? file(params.heme_db) : error("Error: If using --annotate, you must supply prebuilt databases. HEME database file not found at ${params.heme_db}")
    // }

    if (params.use_sulfur) {
        ch_sulfur_db = file(params.sulfur_db).exists() ? file(params.sulfur_db) : error("Error: If using --annotate, you must supply prebuilt databases. SULURR database file not found at ${params.sulfur_db}")
    }

    if (params.use_uniref) {
        ch_uniref_db = file(params.uniref_db).exists() ? file(params.uniref_db) : error("Error: If using --annotate, you must supply prebuilt databases. UNIREF database file not found at ${params.uniref_db}")
        index_mmseqs = true
    }

    if (params.use_methyl) {
        ch_methyl_db = file(params.methyl_db).exists() ? file(params.methyl_db) : error("Error: If using --annotate, you must supply prebuilt databases. METHYL database file not found at ${params.methyl_db}")
        index_mmseqs = true
    }

    if (params.use_fegenie) {
        ch_fegenie_db = file(params.fegenie_db).exists() ? file(params.fegenie_db) : error("Error: If using --annotate, you must supply prebuilt databases. FEGENIE database file not found at ${params.fegenie_db}")
    }

    if (params.use_canthyd) {
        ch_canthyd_hmm_db = file(params.canthyd_hmm_db).exists() ? file(params.canthyd_hmm_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD HMM database file not found at ${params.canthyd_hmm_db}")
        ch_canthyd_mmseqs_db = file(params.canthyd_mmseqs_db).exists() ? file(params.canthyd_mmseqs_db) : error("Error: If using --annotate, you must supply prebuilt databases. CANT_HYD MMseqs database file not found at ${params.canthyd_mmseqs_db}")
        index_mmseqs = true
        ch_canthyd_mmseqs_list = file(params.canthyd_mmseqs_list)
    }

    if (params.use_vogdb) {
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
    ch_methyl_db
    ch_fegenie_db
    ch_canthyd_hmm_db
    ch_canthyd_mmseqs_db
    ch_canthyd_mmseqs_list
    ch_vogdb_db
    ch_viral_db
    index_mmseqs
}