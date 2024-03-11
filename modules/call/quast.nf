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
    run_quast_with_conda(fasta_file_paths, 'quast_results', ${params.threads}, conda_env_path, conda_env_name)

    # Read the single QUAST report generated for all samples
    quast_report_path = 'quast_results/report.tsv'
    if os.path.exists(quast_report_path):
        report_df = pd.read_csv(quast_report_path, sep='\\t')

        # Assuming the QUAST report includes columns for each sample
        # Initialize a list to collect data for the combined report
        collected_data = []

        # Loop through all GFF files, count genes
        for gff_file in gff_file_paths:
            # Extract sample name
            sample_name = os.path.basename(gff_file).split('_')[0]
            num_genes = count_genes_in_gff(gff_file)
            
            # Extract metrics from the QUAST report for this sample
            # Replace the following lines with the actual logic to extract data for sample_name
            no_contigs = 'NA' # Placeholder, replace with actual extraction logic
            largest_contig = 'NA' # Placeholder, replace with actual extraction logic
            N50 = 'NA' # Placeholder, replace with actual extraction logic

            collected_data.append({
                'sample': sample_name,
                'no. contigs': no_contigs,
                'largest contig': largest_contig,
                'N50': N50,
                'no. pred. genes': num_genes
            })

        # Create a DataFrame from the collected data and save it as a TSV file
        combined_df = pd.DataFrame(collected_data)
        combined_df.to_csv('collected_quast.tsv', sep='\\t', index=False)
    else:
        print("QUAST report not found.")
    """
}