process QUAST {

    errorStrategy 'finish'

    input:
    path (collected_fasta_gff)

    output:
    path("quast_results/report.tsv"), emit: quast_tsv
    path("quast_results/icarus.html")
    path("quast_results/icarus_viewers/")   
    path("quast_results/report.html")
    path("quast_results/report.pdf"), optional: true
    path("collected_quast.tsv"), emit: quast_collected_out

    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    import subprocess
    from glob import glob

    # Function to run QUAST
    def run_quast(fasta_files, output_dir, threads):
        quast_cmd = f'quast.py -o {output_dir} --threads {threads} ' + ' '.join(fasta_files)
        subprocess.run(quast_cmd, shell=True, check=True, executable='/bin/bash')

    # Function to count predicted genes in a GFF file
    def count_genes_in_gff(gff_file):
        with open(gff_file, 'r') as file:
            return sum(1 for line in file if '\\tCDS\\t' in line)

    # Find all FASTA and GFF files
    fasta_file_paths = glob('*.fa')
    gff_file_paths = glob('*.gff')

    # Run QUAST on all FASTA files
    threads = ${task.cpus}
    run_quast(fasta_file_paths, 'quast_results', threads)

    # Read the QUAST report generated for all samples
    quast_report_path = 'quast_results/report.tsv'
    report_df = pd.read_csv(quast_report_path, sep='\t', index_col='Assembly')

    # Dynamically identify sample names based on FASTA filenames and match with QUAST report
    collected_data = []
    for gff_file in gff_file_paths:
        base_name = os.path.splitext(os.path.basename(gff_file))[0]
        sample_name = base_name.split('_called_genes')[0]
        num_genes = count_genes_in_gff(gff_file)

        # Extract metrics from the QUAST report for this sample
        for column in report_df.columns:
            if sample_name in column:
                metrics = {
                    'sample': sample_name,
                    'assembly length': report_df.loc['Total length', column],
                    'no. contigs': report_df.loc['# contigs', column],
                    'largest contig': report_df.loc['Largest contig', column],
                    'N50': report_df.loc['N50', column],
                    'GC (%)': report_df.loc['GC (%)', column],
                    'no. pred. genes': num_genes
                }
                collected_data.append(metrics)
                break

    # Convert the list of dictionaries into a DataFrame
    combined_df = pd.DataFrame(collected_data)

    # Save the DataFrame to a TSV file
    combined_df.to_csv('collected_quast.tsv', sep='\t', index=False)
    """
}
