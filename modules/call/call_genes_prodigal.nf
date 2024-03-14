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

    env {
        TMPDIR = './tmp/'
    }

    script:
    """
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

        for ext in fna faa gff; do
            if [ -s "${sample}_called_genes_needs_renaming.\$ext" ]; then
                awk '/^>/ { 
                    split($0, parts, "_"); 
                    gene_number = sprintf("%06d", parts[length(parts)]); 
                    parts[length(parts)] = gene_number; 
                    $0 = parts[1]; 
                    for (i = 2; i in parts; i++) { 
                        $0 = $0 "_" parts[i]; 
                    } 
                } { print }' "${sample}_called_genes_needs_renaming.\$ext" > "${sample}_called_genes.\$ext"
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




/*
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
        -o "${sample}_called_genes_needs_renamaing.gff" \\
        -p ${params.prodigal_mode} \\
        -g ${params.prodigal_trans_table} \\
        -f gff \\
        -d "${sample}_called_genes_needs_renamaing.fna" \\
        -a "${sample}_called_genes_needs_renamaing.faa"

        if [ ! -s "${sample}_called_genes_needs_renamaing.gff" ]; then
            echo "Sample ${sample} - Warning: called genes GFF file is empty or does not exist."
            # Consider managing absent output files for downstream processes here
        else
            # implement renaming here 
            # rename files to ${sample}_called_genes.<fna|gff|faa>

            python ${ch_called_genes_loc_script} "${sample}_called_genes.gff" > "${sample}_called_genes_table.tsv"
        fi
    fi
*/