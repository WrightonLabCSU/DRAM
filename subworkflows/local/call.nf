//
// Subworkflow with functionality specific to the BortonWrightonLabs/dram pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALL_GENES as CALL_GENES_SMALL                } from "../../modules/local/call/call_genes_prodigal.nf"
include { CALL_GENES as CALL_GENES_MEDIUM               } from "../../modules/local/call/call_genes_prodigal.nf"
include { CALL_GENES as CALL_GENES_LARGE                } from "../../modules/local/call/call_genes_prodigal.nf"
include { QUAST                                         } from "../../modules/local/call/quast.nf"
include {resourceBytes; resourceClass                   } from './utils_resource_classes.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL PRODIGAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CALL {
    take:
    ch_fasta  // channel: [ val(input_fasta name), path(fasta) ]

    main:

    // Route each input to a process instance whose initial resources are
    // uniform, allowing supported executors to safely construct job arrays.
    ch_fasta_by_resource = ch_fasta
        .map { name, fasta -> tuple(resourceClass(resourceBytes(fasta)), name, fasta) }
        .branch { resource_class, name, fasta ->
            small: resource_class == 'small'
            medium: resource_class == 'medium'
            large: resource_class == 'large'
        }

    CALL_GENES_SMALL(ch_fasta_by_resource.small)
    CALL_GENES_MEDIUM(ch_fasta_by_resource.medium)
    CALL_GENES_LARGE(ch_fasta_by_resource.large)

    ch_called_genes = CALL_GENES_SMALL.out.prodigal_fna
        .mix(CALL_GENES_MEDIUM.out.prodigal_fna, CALL_GENES_LARGE.out.prodigal_fna)
    ch_called_proteins = CALL_GENES_SMALL.out.prodigal_faa
        .mix(CALL_GENES_MEDIUM.out.prodigal_faa, CALL_GENES_LARGE.out.prodigal_faa)
    ch_gene_locs = CALL_GENES_SMALL.out.prodigal_locs_tsv
        .mix(CALL_GENES_MEDIUM.out.prodigal_locs_tsv, CALL_GENES_LARGE.out.prodigal_locs_tsv)
    ch_gene_gff = CALL_GENES_SMALL.out.prodigal_gff
        .mix(CALL_GENES_MEDIUM.out.prodigal_gff, CALL_GENES_LARGE.out.prodigal_gff)
    ch_filtered_fasta = CALL_GENES_SMALL.out.prodigal_filtered_fasta
        .mix(CALL_GENES_MEDIUM.out.prodigal_filtered_fasta, CALL_GENES_LARGE.out.prodigal_filtered_fasta)

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
        .mix( ch_filtered_fasta.map { _name, file -> file}, ch_gene_gff.map {_name, file -> file} )
        .collect()
        .set { ch_collected_fasta }

    // Run QUAST on individual FASTA file combined with respective GFF
    ch_quast_input = ch_collected_fasta
        .map { files -> tuple(resourceClass(resourceBytes(files)), files) }
    QUAST(ch_quast_input)
    ch_quast_stats = QUAST.out.quast_collected_out

    emit:
    ch_quast_stats
    ch_gene_locs  // channel: [ val(input_fasta name), path(gene_locs_tsv) ]
    ch_called_genes // channel: [ val(input_fasta name), path(called_genes_file.fna) ]
    ch_called_proteins  // channel: [ val(input_fasta name), path(called_proteins_file.faa) ]
    ch_collected_faa
    ch_collected_fna
    ch_collected_fasta
    ch_gene_gff
    ch_filtered_fasta
}
