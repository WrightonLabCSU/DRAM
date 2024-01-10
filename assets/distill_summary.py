import argparse
import pandas as pd

def distill_summary(combined_annotations, genome_summary_form, output_file):
    # Read the input files using pandas
    combined_data = pd.read_csv(combined_annotations, sep='\t')
    genome_summary_data = pd.read_csv(genome_summary_form, sep='\t')

    # Select relevant columns from genome_summary_data for the output
    output_columns = [col for col in genome_summary_data.columns if col != 'gene_id']

    # Merge data based on 'gene_id' columns
    merged_data = pd.merge(combined_data, genome_summary_data, left_on='gene_id', right_on='gene_id', how='inner')

    # Write the distilled summary to the output file
    merged_data[output_columns].to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    # Call the distill_summary function with the provided arguments
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
