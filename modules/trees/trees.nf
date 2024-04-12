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
    #mkdir -p protein_fastas
    #for faa in \$(ls \${ch_collected_proteins}); do
    #    ln -s \${faa} protein_fastas/
    #done

    KO_LIST=${tree_option == 'nar_nxr' ? nar_nxr_ko_list : amoa_pmoa_ko_list}
    python ${ch_trees_scripts}/parse_annotations.py ${annotations_sqlite3} \${KO_LIST} "extracted_query_ids.txt"

    # Loop through each line in the output file, extract the corresponding sequence
    mkdir -p extracted_sequences
    while read sample query_id; do
        seqtk subseq protein_fastas/\${sample}_called_genes.faa <(echo \${query_id}) > extracted_sequences/\${sample}_\${query_id}.fasta
    done < extracted_query_ids.txt

    cat extracted_sequences/*.fasta > combined_extracted_sequences.fasta
    #pplacer -c ${tree_data_files}/${tree_option}/${tree_option}.refpkg combined_extracted_sequences.fasta


    """
}