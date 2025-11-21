include { RENAME_FASTA           } from "../../modules/local/rename/rename_fasta.nf"
include { CALL                   } from "../../subworkflows/local/call.nf"
include { QC                     } from "../../subworkflows/local/qc.nf"
include { DB_SEARCH              } from "../../subworkflows/local/db_search.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO ANNOTATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    ch_fasta  // channel: [ val(input_fasta name), path(fasta) ]
    default_sheet // Path to dummy sheet
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
    use_vog

    main:
    n_fastas = 0
    ch_rrna_collected = default_sheet
    ch_trna_collected = default_sheet
    ch_combined_annotations = default_sheet

    if (params.rename || call) {
        fasta_name = ch_fasta.map { it[0] }
        fasta_files = ch_fasta.map { it[1] }

        n_fastas = file("$params.input_fasta/${params.fasta_fmt}").size()
    }

    if( params.rename ) {
        // We need to use collect so that we pass all the fasta files to the rename process at once
        // Otherwise, it will try to rename each fasta file one at a time
        // Which since rename is so fast, will clog up job queues
        // so it is faster to rename all at once
        RENAME_FASTA( fasta_name.toList(), fasta_files.toList() )
        // we use flatten here to turn a list back into a channel
        renamed_fasta_paths = RENAME_FASTA.out.renamed_fasta_paths.flatten()
        // we need to recreate the fasta channel with the renamed fasta files
        ch_fasta = renamed_fasta_paths.map {
            fasta_name = it.getBaseName()
            tuple(fasta_name, it)
        }
    }

    ch_quast_stats = default_sheet
    ch_gene_locs = default_sheet
    ch_called_proteins = default_sheet
    ch_collected_fna = default_sheet

    if (call){
        CALL( ch_fasta )
        ch_quast_stats = CALL.out.ch_quast_stats
        ch_gene_locs = CALL.out.ch_gene_locs
        ch_called_proteins = CALL.out.ch_called_proteins
        ch_collected_fna = CALL.out.ch_collected_fna

    }

 
    if (params.annotate){
        DB_SEARCH( 
            ch_gene_locs, 
            ch_called_proteins, 
            default_sheet, 
            n_fastas, 
            call,
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
            use_vog
            )
        ch_combined_annotations = DB_SEARCH.out.ch_combined_annotations
    }

    if (params.qc){
        QC( ch_fasta, default_sheet, ch_combined_annotations, ch_collected_fna, call )
        ch_rrna_collected = QC.out.ch_rrna_collected
        ch_trna_collected = QC.out.ch_trna_collected
        ch_combined_annotations = QC.out.ch_final_annots
    }

    emit:
    ch_rrna_collected
    ch_trna_collected
    ch_combined_annotations
    ch_quast_stats

}