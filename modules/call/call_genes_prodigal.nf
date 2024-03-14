process CALL_GENES {
    // Your process configuration remains the same...

    script:
    """
    mkdir -p tmp
    export TMPDIR='./tmp/'

    /opt/bbmap/reformat.sh \\
    in=${fasta} \\
    out="${sample}_${params.min_contig_len}.fa" \\
    minlength=${params.min_contig_len}

    if [ ! -s "${sample}_${params.min_contig_len}.fa" ]; then
        echo "Sample ${sample} - Error: no contigs after filtering for minimum contig length (${params.min_contig_len})"
        exit 1
    else
        prodigal \\
        -i "${sample}_${params.min_contig_len}.fa" \\
        -o "${sample}_called_genes_needs_renaming.gff" \\
        -d "${sample}_called_genes_needs_renaming.fna" \\
        -a "${sample}_called_genes_needs_renaming.faa" \\
        -p ${params.prodigal_mode} \\
        -g ${params.prodigal_trans_table} \\
        -f gff

        gene_counter=1

        # Process .fna and .faa
        for ext in fna faa; do
            awk -v prefix=">${sample}_" '/^>/{ 
                print prefix sprintf("%06d", gene_counter); 
                gene_counter++;
                next; 
            }
            { print }' "${sample}_called_genes_needs_renaming.\$ext" > "${sample}_called_genes.\$ext"
        done

        # Reset counter for .gff processing
        gene_counter=1

        # Process .gff
        awk -v prefix="${sample}_" 'BEGIN{FS=OFS="\t"}
            !/^#/ && \$3 == "CDS" {
                gsub(/ID=[^;]+/, "ID=" prefix sprintf("%06d", gene_counter), \$9);
                gene_counter++;
            }
            {print}' "${sample}_called_genes_needs_renaming.gff" > "${sample}_called_genes.gff"

        # Now, process the .gff file to generate .tsv using the Python script
        if [ -s "${sample}_called_genes.gff" ]; then
            python ${ch_called_genes_loc_script} "${sample}_called_genes.gff" > "${sample}_called_genes_table.tsv"
        else
            echo "Sample ${sample} - Warning: called genes GFF file is empty or does not exist."
        fi
    fi
    """
}
