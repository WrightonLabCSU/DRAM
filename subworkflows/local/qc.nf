include { COLLECT_RNA            } from "../../subworkflows/local/collect_rna.nf"
include { ADD_TAXA               } from "../../modules/local/add_and_combine/add_taxa.nf"
include { ADD_BIN_QUALITY        } from "../../modules/local/add_and_combine/add_bin_quality.nf"
include { GENERATE_GFF_GENBANK   } from "../../modules/local/add_and_combine/generate_gff_genbank.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QC {
    take:
    ch_fasta  // channel: [ val(input_fasta name), path(fasta) ]
    default_sheet // Path to dummy sheet
    ch_combined_annotations  // channel: [ path(combined_annotations_out) ]
    ch_collected_fna
    call                     // boolean: whether gene calling flag is set

    main:

    COLLECT_RNA( ch_fasta, default_sheet, call )
    ch_rrna_collected = COLLECT_RNA.out.ch_rrna_collected
    ch_trna_collected = COLLECT_RNA.out.ch_trna_collected


    // Add Bin Quality to annotations
    if( params.bin_quality ){
        ch_bin_quality = file(params.bin_quality)
        ADD_BIN_QUALITY( ch_combined_annotations, ch_bin_quality )
        ch_updated_annots = ADD_BIN_QUALITY.out.annots_bin_quality_out
    }
    else{
        ch_updated_annots = ch_combined_annotations
    }

    // Add Taxonomy to annotations
    if( params.taxa ){
        ch_taxa = file(params.taxa)
        ADD_TAXA( ch_updated_annots, ch_taxa )
        ch_updated_taxa_annots = ADD_TAXA.out.annots_taxa_out
    }
    else{
        ch_updated_taxa_annots = ch_combined_annotations
    }

    ch_final_annots = ch_updated_taxa_annots


    if( params.generate_gff || params.generate_gbk ){
        if (!call) {
            ch_called_genes = Channel
                .fromPath(file(params.input_genes) / params.genes_fna_fmt, checkIfExists: true)
                .ifEmpty { exit 1, "If you specify --generate_gff or --generate_gbk without --call, you must provide a fasta file of called genes using --input_genes and --genes_fna_fmt,. Cannot find any called gene fasta files matching: ${params.input_genes} and ${params.genes_fna_fmt}\nNB: Path needs to follow pattern: path/to/directory/" }
                .map {
                    input_fastaName = it.getBaseName()
                    tuple(input_fastaName, it)
                }
            // Collect all individual fasta to pass to quast
            Channel.empty()
                .mix( ch_called_genes  )
                .collect()
                .set { ch_collected_fna }
        }
    
        GENERATE_GFF_GENBANK( ch_collected_fna, params.database_list, ch_final_annots )
    }




    emit:
    ch_final_annots
    ch_rrna_collected
    ch_trna_collected
}