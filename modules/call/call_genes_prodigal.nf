process CALL_GENES {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)
    path(ch_called_genes_loc_script)

    output:
    tuple val(sample), path("${sample}_called_genes.fna"), emit: prodigal_fna, optional: true
    tuple val(sample), path("${sample}_called_genes.faa"), emit: prodigal_faa, optional: true
    tuple val(sample), path("${sample}_called_genes_table.tsv"), emit: prodigal_locs_tsv, optional: true
    path("${sample}_${params.min_contig_len}.fa"), emit: prodigal_filtered_fasta, optional: true
    path("${sample}_called_genes.gff"), emit: prodigal_gff, optional: true

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

        gene_counter=1 # Initialize the gene counter correctly

        for ext in fna faa gff; do
            if [[ "$ext" == "gff" ]]; then
                awk -v prefix="${sample}_" 'BEGIN{FS=OFS="\t"; gene_counter=1}
                    $0 !~ /^#/ {
                        ...
                        # Processing lines, incrementing gene_counter
                        gene_counter++;
                    }
                    {print}' "${sample}_called_genes_needs_renaming.$ext" > "${sample}_called_genes.$ext"
            else
                # For .fna and .faa
                awk -v prefix=">${sample}_" '/^>/{ 
                    ...
                    print prefix sprintf("%06d", gene_counter); # Correct use before increment
                    gene_counter++;
                    next; 
                }
                { print }' "${sample}_called_genes_needs_renaming.$ext" > "${sample}_called_genes.$ext"
            fi
        done

        if [ -s "${sample}_called_genes.gff" ]; then
            python ${ch_called_genes_loc_script} "${sample}_called_genes.gff" > "${sample}_called_genes_table.tsv"
        else
            echo "Sample ${sample} - Warning: called genes GFF file is empty or does not exist."
        fi
    fi
    """
}
