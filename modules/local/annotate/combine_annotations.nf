process COMBINE_ANNOTATIONS {
    label 'process_low'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_hmmer_mmseqs2_pruned:fc59f737a5e0566a"

    input:
    val all_annotations
    val all_genes
    // tuple val( input_fasta ), path( fasta )

    output:
    path "raw-annotations.tsv", emit: combined_annotations_out

    script:
    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    # Create a log directory if it doesn't exist
    mkdir logs

    # Define the log file path
    touch logs/combine_annotations.log
    log_file="logs/combine_annotations.log"

    combine_annotations.py --annotations ${all_annotations} --threads ${params.threads} --output "raw-annotations.tsv" >> \$log_file 2>&1 --genes-faa ${all_genes}

    """
}
