import argparse
import pandas as pd

def distill_summary(combined_annotations_path, genome_summary_form_path, output_path):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_path, sep='\t')

    # Initialize the output DataFrame
    distill_summary_df = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    # Iterate through gene_id values in genome_summary_form
    for gene_id in genome_summary_form['gene_id']:
        # Search for a match in combined_annotations columns ending in "_id"
        match_columns = [col for col in combined_annotations.columns if col.endswith('_id') and col != 'query_id']

        for column in match_columns:
            # Check for a match
            match_rows = combined_annotations[combined_annotations[column] == gene_id]

            if not match_rows.empty:
                # Add matching information to the output DataFrame
                distill_summary_df = distill_summary_df.append({
                    'gene_id': gene_id,
                    'query_id': match_rows['query_id'].values[0],
                    'sample': match_rows['sample'].values[0],
                }, ignore_index=True)

    # Write the output to a TSV file
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    # Call the distill_summary function with provided arguments
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
