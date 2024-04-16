process TREES {

    errorStrategy 'finish'

    input:
    path( ch_combined_annotations, stageAs: "raw-annotations.tsv" )
    val( tree_option )
    path( ch_collected_proteins )
    path( tree_data_files )
    path( ch_trees_scripts )
    val( nar_nxr_ko_list )
    val( amoa_pmoa_ko_list )

    output:
    path("updated-annotations.tsv"), emit: trees_out, optional: true

    script:
    """

    ln -s ${tree_data_files}/* .
    ln -s ${ch_trees_scripts}/*.py .

    KO_LIST="${tree_option == 'nar_nxr' ? nar_nxr_ko_list : amoa_pmoa_ko_list}"
    python parse_annotations.py raw-annotations.tsv \${KO_LIST} "extracted_query_ids.txt"

    if [ -s "extracted_query_ids.txt" ]; then
        # The file is not empty, proceed with processing

        # Loop through each line in the output file, extract the corresponding sequence
        mkdir -p extracted_sequences
        while IFS=\$'\t' read -r sample query_id; do
            seqtk subseq \${sample}_called_genes.faa <(echo \${query_id}) > extracted_sequences/\${sample}_\${query_id}.fasta
        done < extracted_query_ids.txt

        # Combine all sequences into one file
        cat extracted_sequences/*.fasta > combined_extracted_sequences.fasta

        # Align sequences to the reference alignment
        mafft --thread ${params.threads} --add combined_extracted_sequences.fasta --reorder trees/${tree_option}/${tree_option}.refpkg/nar_nxr_seqs_for_tree_aligned.faa > aligned_sequences.fasta

        cat extracted_sequences/*.fasta > combined_extracted_sequences.fasta
        # Uncomment the following line to run pplacer if the rest of the script works fine
        pplacer -j ${task.cpus} -c trees/${tree_option}/${tree_option}.refpkg aligned_sequences.fasta

        python update_annots_trees.py aligned_sequences.jplace raw-annotations.tsv ${tree_option}/${tree_option}-tree-mapping.tsv updated-annotations.tsv
    else
        echo "No gene IDs of interest found, skipping sequence extraction and analysis."
    fi

    """

}
