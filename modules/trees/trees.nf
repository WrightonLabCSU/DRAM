process TREES {

    errorStrategy 'finish'

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    path( annotations_sqlite3 )
    val( tree_option )
    path( ch_collected_proteins )
    path( tree_data_files )
    path( ch_trees_scripts )
    val( nar_nxr_ko_list )
    val( amoa_pmoa_ko_list )



    output:
    path("something.tre"), emit: trees_out, optional: true

    script:
    """

    KO_LIST="${tree_option == 'nar_nxr' ? nar_nxr_ko_list : amoa_pmoa_ko_list}"
    python ${ch_trees_scripts}/parse_annotations.py ${annotations_sqlite3} \${KO_LIST} "extracted_query_ids.txt"

    # Loop through each line in the output file, extract the corresponding sequence
    mkdir -p extracted_sequences
    while IFS=$'\t' read -r sample query_id; do
        seqtk subseq \${sample}_called_genes.faa <(echo \${query_id}) > extracted_sequences/\${sample}_\${query_id}.fasta
    done < extracted_query_ids.txt

    cat extracted_sequences/*.fasta > combined_extracted_sequences.fasta
    # Uncomment the following line to run pplacer if the rest of the script works fine
    # pplacer -c \${tree_data_files}/\${tree_option}/\${tree_option}.refpkg combined_extracted_sequences.fasta
    """

}