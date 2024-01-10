import argparse
import pandas as pd

def distill_summary(combined_annotations_file, genome_summary_form_file, output_file):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_file, sep='\t')

    # Print gene_ids in genome_summary_form for debugging
    print("gene_ids in genome_summary_form:")
    print(genome_summary_form['gene_id'])

    # Print columns in combined_annotations for debugging
    print("Columns in combined_annotations:")
    print(combined_annotations.columns)

    # Merge DataFrames based on gene_id and query_id
    merged_data = pd.merge(genome_summary_form, combined_annotations, left_on='gene_id', right_on='query_id', how='inner')

    # Print gene_ids after merging for debugging
    print("gene_ids after merging:")
    print(merged_data['gene_id'].unique())

    # Extract relevant columns
    distilled_summary = merged_data[['gene_id', 'query_id', 'sample']]

    # Print distilled summary for debugging
    print("Distilled Summary:")
    print(distilled_summary)

    # Write the distilled summary to the output file
    distilled_summary.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.genome_summary_form, args.output)
