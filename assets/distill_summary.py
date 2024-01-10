import argparse
import pandas as pd

def distill_summary(combined_annotations_path, genome_summary_form_path, output_path):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_path, sep='\t')

    # Initialize the output DataFrame
    distill_summary_df = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    # Iterate through gene_ids in genome_summary_form
    for gene_id in genome_summary_form['gene_id']:
        # Find matching rows in combined_annotations based on gene_id
        matching_rows = combined_annotations[combined_annotations[gene_id + '_id'] == gene_id]

        # Add information to distill_summary_df if matches are found
        if not matching_rows.empty:
            distill_summary_df = distill_summary_df.append({
                'gene_id': gene_id,
                'query_id': matching_rows['query_id'].iloc[0],
                'sample': matching_rows['sample'].iloc[0]
            }, ignore_index=True)

    # Write the distilled summary to the output file
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
