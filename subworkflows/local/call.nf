//
// Subworkflow with functionality specific to the WrightonLabCSU/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALL_GENES                                    } from "../../modules/local/call/call_genes_prodigal.nf"
include { QUAST                                         } from "../../modules/local/call/quast.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL PRODIGAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CALL {
    take:
    ch_fasta  // channel: [ val(input_fasta name), path(fasta) ]

    main:

    // Call genes using Prodigal on the input fasta file(s) 1-by-1
    CALL_GENES ( ch_fasta )
    ch_called_genes = CALL_GENES.out.prodigal_fna
    ch_called_proteins = CALL_GENES.out.prodigal_faa
    ch_gene_locs = CALL_GENES.out.prodigal_locs_tsv
    ch_gene_gff = CALL_GENES.out.prodigal_gff
    ch_filtered_fasta = CALL_GENES.out.prodigal_filtered_fasta

    // Collect all individual fasta to pass to quast
    ch_called_proteins
        .map { tuple -> tuple[1] }  // Extract only the file path from each tuple
        .collect()                  // Collect all paths into a list
        .set { ch_collected_faa }   // Set the resulting list to ch_collected_faa

    // Collect all individual fasta to pass to quast
    channel.empty()
        .mix( ch_called_genes  )
        .collect()
        .set { ch_collected_fna }

    // Collect all individual fasta to pass to quast
    channel.empty()
        .mix( ch_filtered_fasta, ch_gene_gff  )
        .collect()
        .set { ch_collected_fasta }

    // QUAST is per-input_fasta; for DRAM-v viral mode the unit of analysis is
    // the scaffold (per-vMAG), and the catalog-level QUAST stats would collapse
    // all contigs into one row. Skip the process and substitute the dummy sheet.
    if (!params.use_dramv) {
        QUAST( ch_collected_fasta )
        ch_quast_stats = QUAST.out.quast_collected_out
    } else {
        ch_quast_stats = Channel.value(file(params.distill_dummy_sheet))
    }

    emit:
    ch_quast_stats
    ch_gene_locs  // channel: [ val(input_fasta name), path(gene_locs_tsv) ]
    ch_called_genes // channel: [ val(input_fasta name), path(called_genes_file.fna) ]
    ch_called_proteins  // channel: [ val(input_fasta name), path(called_proteins_file.faa) ]
    ch_collected_faa  
    ch_collected_fna
    ch_collected_fasta
}
