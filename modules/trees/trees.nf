process TREES {

    errorStrategy 'finish'

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    val( tree_option )
    path( ch_collected_proteins, stageAs: "protein_fastas/*" )
    path( tree_data_files )
    path( ch_trees_scripts )



    output:
    path("something.tre"), emit: trees_out, optional: true

    script:
    """
    # Call the Python script to parse annotations and extract query IDs
    python ${ch_trees_scripts}/parse_annotations.py ${ch_combined_annotations} ${tree_option} "extracted_query_ids.txt"

    # Extract sequences based on query IDs
    seqtk subseq ${ch_collected_proteins} extracted_query_ids.txt > extracted_sequences.fasta

    # Run pplacer (assuming pplacer command and options need to be adjusted based on your setup)
    pplacer -c ${tree_data_files}/${tree_option}/${tree_option}.refpkg extracted_sequences.fasta

    # Assuming pplacer generates a tree, adjust commands and filenames as necessary


    """

}