import pandas as pd
import argparse

def distill_summary(combined_annotations, target_id_counts, genome_summary_form, output_file):
    # Load dataframes from TSV files
    combined_df = pd.read_csv(combined_annotations, sep='\t')
    target_id_counts_df = pd.read_csv(target_id_counts, sep='\t')
    genome_summary_form_df = pd.read_csv(genome_summary_form, sep='\t')

    # Extract columns ending with '_id' from combined_annotations
    id_columns = [col for col in combined_df.columns if col.endswith('_id')]
    gene_id_column = 'gene_id'
    
    # Create a new dataframe with required columns
    distill_df = combined_df[['query_id', 'sample'] + id_columns].copy()

    # Merge with target_id_counts to get sample columns
    distill_df = pd.merge(distill_df, target_id_counts_df, how='left', left_on='sample', right_on='target_id')

    # Merge with genome_summary_form using gene_id
    distill_df = pd.merge(distill_df, genome_summary_form_df, how='left', left_on=gene_id_column, right_on='gene_id')

    # Drop unnecessary columns and write to output file
    distill_df.drop(columns=['target_id', gene_id_column], inplace=True)
    distill_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Distill summary script')
    parser.add_argument('--combined_annotations', required=True, help='Path to combined_annotations TSV file')
    parser.add_argument('--target_id_counts', required=True, help='Path to target_id_counts TSV file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to genome_summary_form TSV file')
    parser.add_argument('--output', required=True, help='Path to output distill_summary.tsv file')

    args = parser.parse_args()

    distill_summary(args.combined_annotations, args.target_id_counts, args.genome_summary_form, args.output)
