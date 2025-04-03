process CALL_GENES {
    label 'process_low'

    errorStrategy 'finish'

    tag { input_fasta }

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_scikit-bio_hmmer_pruned:ef64c488c99048d6"

    input:
    tuple val( input_fasta ), path( fasta )

    output:
    tuple val( input_fasta ), path( "${input_fasta}_called_genes.fna" ), emit: prodigal_fna, optional: true
    tuple val( input_fasta ), path( "${input_fasta}_called_genes.faa" ), emit: prodigal_faa, optional: true
    tuple val( input_fasta ), path( "${input_fasta}_called_genes_table.tsv" ), emit: prodigal_locs_tsv, optional: true
    path( "${input_fasta}_${params.min_contig_len}.fa" ), emit: prodigal_filtered_fasta, optional: true
    path( "${input_fasta}_called_genes.gff" ), emit: prodigal_gff, optional: true


    script:

    """

    reformat.sh \\
    in=${fasta} \\
    out="${input_fasta}_${params.min_contig_len}.fa" \\
    minlength=${params.min_contig_len}

    if [ ! -s "${input_fasta}_${params.min_contig_len}.fa" ]; then
        echo "Sample ${input_fasta} - Error: no contigs after filtering for minimum contig length (${params.min_contig_len})"
        exit 1
    else
        prodigal \\
        -i "${input_fasta}_${params.min_contig_len}.fa" \\
        -o "${input_fasta}_called_genes.gff" \\
        -p ${params.prodigal_mode} \\
        -g ${params.prodigal_trans_table} \\
        -f gff \\
        -d "${input_fasta}_called_genes.fna" \\
        -a "${input_fasta}_called_genes.faa"

        if [ ! -s "${input_fasta}_called_genes.gff" ]; then
            echo "Sample ${input_fasta} - Warning: called genes GFF file is empty or does not exist."
            # Consider managing absent output files for downstream processes here
        else
            generate_gene_loc_tsv.py "${input_fasta}_called_genes.gff" > "${input_fasta}_called_genes_table.tsv"
        fi
    fi
    """
}
