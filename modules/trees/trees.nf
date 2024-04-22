process TREES {

    errorStrategy 'finish'

    input:
    path( ch_combined_annotations, stageAs: "initial-annotations.tsv" )
    val( trees_list )
    path( ch_collected_proteins )
    path( tree_data_files )
    path( ch_trees_scripts )
    file( ch_add_trees )

    output:
    path("updated-annotations.tsv"), emit: updated_annotations, optional: true

    script:
    """        
    ln -s ${tree_data_files}/* .
    ln -s ${ch_trees_scripts}/*.py .

    cp initial-annotations.tsv current-annotations.tsv

    # Symlink additional tree directories if provided
    if [[ "${params.add_trees}" != "" ]]; then
        ln -s ${ch_add_trees}/* ${tree_data_files}/
    fi

    # Append additional tree directories to trees_list
    if [ -d "${params.add_trees}" ]; then
        for dir in \$(ls ${params.add_trees}); do
            trees_list+=";\$dir"
        done
    fi

    # Split trees_list into an array
    IFS=';' read -ra TREE_OPTIONS <<< "${trees_list}"
    
    for tree_option in "\${TREE_OPTIONS[@]}"; do
        echo "Processing tree: \${tree_option}"

        # Determine the search terms file for the current tree
        KO_LIST="\${tree_option}/\${tree_option}_search_terms.txt"
        python parse_annotations.py current-annotations.tsv \${KO_LIST} "extracted_query_ids.txt"
        
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
            mafft --thread ${params.threads} --add combined_extracted_sequences.fasta --reorder trees/\${tree_option}/\${tree_option}.refpkg/nar_nxr_seqs_for_tree_aligned.faa > aligned_sequences.fasta
            
            # Run pplacer
            pplacer -j ${task.cpus} -c trees/\${tree_option}/\${tree_option}.refpkg aligned_sequences.fasta
            
            # Update the annotations using the mapping and the placements
            python update_annots_trees.py aligned_sequences.jplace current-annotations.tsv "trees/\${tree_option}/\${tree_option}-tree-mapping.tsv" updated-annotations.tsv

            # Set the updated annotations as the current for the next tree
            mv updated-annotations.tsv current-annotations.tsv
        else
            echo "No gene IDs of interest found for tree \${tree_option}, skipping sequence extraction and analysis."
        fi
    done

    # Finalize the process
    mv current-annotations.tsv updated-annotations.tsv

    """
}