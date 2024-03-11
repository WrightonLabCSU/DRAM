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

    # Dynamically identify sample names based on FASTA filenames
    sample_names = [os.path.splitext(os.path.basename(fasta))[0] for fasta in fasta_file_paths]

    # Loop through all GFF files, count genes, and match with QUAST report
    collected_data = []
    for gff_file in gff_file_paths:
        base_name = os.path.splitext(os.path.basename(gff_file))[0]
        sample_name = base_name.split('_called_genes')[0]
        num_genes = count_genes_in_gff(gff_file)

        # Find the corresponding column in the QUAST report for this sample
        for column in report_df.columns:
            if sample_name in column:
                metrics = report_df[column].to_dict()
                metrics['sample'] = sample_name
                metrics['no. pred. genes'] = num_genes
                collected_data.append(metrics)
                break

    # Create a DataFrame from the collected data
    combined_df = pd.DataFrame(collected_data)

    # If needed, rename or select specific columns from the combined DataFrame
    # Example renaming, adjust according to actual column names found in report_df
    combined_df.rename(columns={
        'Total length': 'total length',
        'Largest contig': 'largest contig',
        'N50': 'N50',
        'no. pred. genes': 'no. pred. genes'
    }, inplace=True)

    # Select and reorder columns based on your specific needs
    desired_columns = ['sample', 'total length', 'largest contig', 'N50', 'no. pred. genes']
    combined_df = combined_df[desired_columns]

    # Save the DataFrame to a TSV file
    combined_df.to_csv('collected_quast.tsv', sep='\t', index=False)
    """
}