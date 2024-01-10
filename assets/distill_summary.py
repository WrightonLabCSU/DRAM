import argparse
import pandas as pd

def distill_summary(combined_annotations_file, genome_summary_file, output_file):
    # Read input files into DataFrames
    combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
    genome_summary = pd.read_csv(genome_summary_file, sep='\t')

    # Initialize the output DataFrame
    distill_summary_df = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    # Iterate through gene_ids in genome_summary_file
    for gene_id in genome_summary['gene_id']:
        # Find matching row in combined_annotations based on _id columns
        match_row = combined_annotations[combined_annotations.filter(like='_id').eq(gene_id).any(axis=1)]

        if not match_row.empty:
            # Extract relevant information and append to distill_summary_df
            query_id = match_row['query_id'].iloc[0]
            sample = match_row['sample'].iloc[0]
            distill_summary_df = distill_summary_df.append({'gene_id': gene_id, 'query_id': query_id, 'sample': sample}, ignore_index=True)

    # Write the result to the output file
    distill_summary_df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    # Call the distill_summary function with the provided arguments
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
