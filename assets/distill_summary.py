import argparse
import pandas as pd

def distill_summary(combined_annotations_file, genome_summary_form_file, output_file):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_file, sep='\t')

    # Extract gene_id, query_id, and sample information
    distilled_summary = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    for index, row in genome_summary_form.iterrows():
        gene_id = row['gene_id']

        # Check if the gene_id is present in combined_annotations
        if f"{gene_id}_id" in combined_annotations.columns:
            matching_data = combined_annotations[['query_id', 'sample', f"{gene_id}_id"]]
            matching_data = matching_data.rename(columns={f"{gene_id}_id": 'query_id'})
            matching_data['gene_id'] = gene_id

            # Append matching data to the distilled_summary dataframe
            distilled_summary = pd.concat([distilled_summary, matching_data], ignore_index=True)

    # Write the distilled summary to the output file
    distilled_summary.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
