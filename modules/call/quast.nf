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

    # Function to run QUAST
    def run_quast(fasta_files, output_dir, threads):
        cmd = ['quast.py', '-o', output_dir, '--threads', str(threads)] + fasta_files
        subprocess.run(cmd, check=True)

    # Function to count predicted genes in a GFF file
    def count_genes_in_gff(gff_file):
        with open(gff_file, 'r') as file:
            return sum(1 for line in file if '\\tgene\\t' in line)

    # Get all fasta and gff file paths
    fasta_file_paths = glob('*.fa')
    gff_file_paths = glob('*.gff')

    # Run QUAST on all collected FASTA files together
    run_quast(fasta_file_paths, 'quast_results', ${params.threads})

    # Initialize a list to collect data for the combined report
    collected_data = []

    # Loop through all GFF files, count genes, and prepare the data for the combined report
    for gff_file in gff_file_paths:
        # Extract sample name by removing '_called_genes.gff'
        sample_name = os.path.basename(gff_file).replace('_called_genes.gff', '')
        num_genes = count_genes_in_gff(gff_file)
        
        # Construct the expected QUAST report path based on the sample name
        quast_report_path = os.path.join('quast_results', f'{sample_name}_report.tsv')
        
        # Check if the report exists before attempting to read it
        if os.path.exists(quast_report_path):
            df = pd.read_csv(quast_report_path, sep='\\t', names=["Metric", "Value"], skiprows=1)
            df.set_index("Metric", inplace=True)
            
            collected_data.append({
                'sample': sample_name,
                'no. contigs': df.at['# contigs', 'Value'],
                'largest contig': df.at['Largest contig', 'Value'],
                'N50': df.at['N50', 'Value'],
                'no. pred. genes': num_genes
            })
        else:
            print(f"Report for {sample_name} not found.")

    # Create a DataFrame from the collected data and save it as a TSV file
    combined_df = pd.DataFrame(collected_data)
    combined_df.to_csv('collected_quast.tsv', sep='\\t', index=False)
    """
}