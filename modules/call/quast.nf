process QUAST {

    errorStrategy 'finish'

    input:
    path ( collected_fasta_gff )

    output:
    path( "quast_results/report.tsv" ), emit: quast_tsv
    path( "${sample}_QUAST/icarus.tsv" )
    path( "${sample}_QUAST/report.html" )
    path( "${sample}_QUAST/report.pdf" )
    path( "collected_quast.tsv" ), emit: quast_collected_out

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    import subprocess
    from glob import glob

    # Function to activate conda environment and run QUAST
    def run_quast_with_conda(fasta_files, output_dir, threads, conda_env_path, conda_env_name):
        activate_env = f'source {conda_env_path}/bin/activate {conda_env_name}'
        quast_cmd = f'quast.py -o {output_dir} --threads {threads} ' + ' '.join(fasta_files)
        cmd = f'{activate_env} && {quast_cmd}'
        subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

    # Function to count predicted genes in a GFF file
    def count_genes_in_gff(gff_file):
        with open(gff_file, 'r') as file:
            return sum(1 for line in file if '\\tgene\\t' in line)

    # Find all FASTA and GFF files in the current directory
    fasta_file_paths = glob('*.fa')
    gff_file_paths = glob('*.gff')

    # Activate conda environment and run QUAST on all FASTA files together
    conda_env_path = '/opt/miniconda'
    conda_env_name = 'support'
    run_quast_with_conda(fasta_file_paths, 'quast_results', '4', conda_env_path, conda_env_name)

    # Read the QUAST report generated for all samples
    quast_report_path = 'quast_results/report.tsv'
    report_df = pd.read_csv(quast_report_path, sep='\t', index_col='Assembly')

    # Loop through all GFF files, count genes, and prepare the data for the combined report
    collected_data = []
    for gff_file in gff_file_paths:
        sample_name = os.path.basename(gff_file).split('_')[0]
        num_genes = count_genes_in_gff(gff_file)
        
        # Extract metrics from the QUAST report for this sample
        metrics = report_df[sample_name + '_2500'].to_dict()
        metrics['sample'] = sample_name
        metrics['no. pred. genes'] = num_genes
        collected_data.append(metrics)

    # Convert the list of dictionaries into a DataFrame
    combined_df = pd.DataFrame(collected_data)

    # If needed, adjust the columns to match your desired output structure
    # For instance, you might only want to keep certain metrics:
    combined_df = combined_df[['sample', 'no. contigs (>= 0 bp)', 'Largest contig', 'N50', 'no. pred. genes']]

    # Rename columns to remove the condition '(>= 0 bp)' for clarity, if desired
    combined_df.rename(columns={
        'no. contigs (>= 0 bp)': 'no. contigs',
        'Largest contig': 'largest contig',
        'N50': 'N50',
        'no. pred. genes': 'no. pred. genes'
    }, inplace=True)

    # Save the DataFrame to a TSV file
    combined_df.to_csv('collected_quast.tsv', sep='\t', index=False)
    """
}