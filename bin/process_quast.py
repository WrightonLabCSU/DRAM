#!/usr/bin/env python
import os
import pandas as pd
from glob import glob
from utils.logger import get_logger
from pathlib import Path

FASTA_COLUMN = os.getenv('FASTA_COLUMN', 'input_fasta')

logger = get_logger(filename=Path(__file__).stem)


# Function to count predicted genes in a GFF file
def count_genes_in_gff(gff_file):
    with open(gff_file, 'r') as file:
        return sum(1 for line in file if '\\tCDS\\t' in line)

# Find all FASTA and GFF files
gff_file_paths = glob('*.gff')

# Read the QUAST report generated for all input_fastas
quast_report_path = 'quast_results/report.tsv'
report_df = pd.read_csv(quast_report_path, sep='\t', index_col='Assembly')

logger.info(f"Number of gff files to process: {len(gff_file_paths)}")

# Dynamically identify input_fasta names based on FASTA filenames and match with QUAST report
collected_data = []
for i, gff_file in enumerate(gff_file_paths):
    try: 
        base_name = os.path.splitext(os.path.basename(gff_file))[0]
        input_fasta_name = base_name.split('_called_genes')[0]
        num_genes = count_genes_in_gff(gff_file)

        # Extract metrics from the QUAST report for this input_fasta
        for column in report_df.columns:
            if input_fasta_name in column:
                metrics = {
                    FASTA_COLUMN: input_fasta_name,
                    'assembly size': report_df.loc['Total length', column],
                    'no. contigs': report_df.loc['# contigs', column],
                    'largest contig': report_df.loc['Largest contig', column],
                    'N50': report_df.loc['N50', column],
                    'L50': report_df.loc['L50', column],
                    'GC (%)': report_df.loc['GC (%)', column],
                    'no. pred. genes': num_genes
                }
                collected_data.append(metrics)
                break
    except Exception as e:
        logger.error(f"Error processing the {i+1}-th GFF file {gff_file}: {str(e)}")
        continue

# Convert the list of dictionaries into a DataFrame
combined_df = pd.DataFrame(collected_data)

# Save the DataFrame to a TSV file
combined_df.to_csv('collected_quast.tsv', sep='\t', index=False)