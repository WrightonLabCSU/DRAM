import argparse
import pandas as pd

def distill_summary(combined_annotations_path, genome_summary_form_path, output_path, add_module_files):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_path, sep='\t')

    # Additional columns from add_moduleX files
    additional_modules = {}
    for i, add_module_file in enumerate(add_module_files, start=1):
        if add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_module_data.columns = [f'{col}_{i}' if col != 'gene_id' else col for col in additional_module_data.columns]
            additional_modules[f'add_module{i}'] = additional_module_data

    # Merge add_moduleX data with genome_summary_form
    for module_name, add_module_data in additional_modules.items():
        # Use explicit suffixes to avoid duplicate columns
        # Exclude 'gene_id' column from add_moduleX files
        cols_to_merge = [col for col in add_module_data.columns if col != 'gene_id']
        genome_summary_form = pd.merge(
            genome_summary_form,
            add_module_data[cols_to_merge],
            left_on='gene_id',
            right_on=f'gene_id_{module_name}',
            how='left'
        )

    # Drop duplicate gene_id columns
    genome_summary_form = genome_summary_form.loc[:,~genome_summary_form.columns.duplicated()]

    # Initialize the output DataFrame
    distill_summary_df = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    # Iterate through gene_id values in genome_summary_form
    for _, row in genome_summary_form.iterrows():
        gene_id = row['gene_id']

        # Check if gene_id is present in combined_annotations columns ending in "_id"
        match_columns = [col for col in combined_annotations.columns if col.endswith('_id') and col != 'query_id']

        if any(combined_annotations[col].eq(gene_id).any() for col in match_columns):
            # If gene_id is present, proceed with adding information to the output DataFrame
            for column in match_columns:
                # Check for matches
                match_rows = combined_annotations[combined_annotations[column] == gene_id]

                if match_rows.empty:
                    print(f"No matches found for gene_id: {gene_id}")  # Debug statement

                for _, match_row in match_rows.iterrows():
                    # Instead of using DataFrame.append, use pandas.concat
                    distill_summary_df = pd.concat([distill_summary_df, pd.DataFrame({
                        'gene_id': gene_id,
                        'query_id': match_row['query_id'],
                        'sample': match_row['sample'],
                    }, index=[0])], ignore_index=True)

                    # Add values from additional columns in genome_summary_form
                    for col in genome_summary_form.columns:
                        if col not in ['gene_id'] + match_columns:
                            distill_summary_df.at[distill_summary_df.index[-1], col] = row[col]

                    # Add values from selected columns in combined_annotations
                    for col in combined_annotations.columns:
                        if col not in ['query_id'] + match_columns:
                            distill_summary_df.at[distill_summary_df.index[-1], col] = match_row[col]

                    # Add values from add_moduleX columns
                    for add_module_col in add_module_data.columns:
                        if add_module_col != 'gene_id':
                            distill_summary_df.at[distill_summary_df.index[-1], add_module_col] = row[add_module_col]

    # Write the output to a TSV file
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')
    parser.add_argument('--add_module1', default='empty', help='Path to add_module1 file')
    parser.add_argument('--add_module2', default='empty', help='Path to add_module2 file')
    parser.add_argument('--add_module3', default='empty', help='Path to add_module3 file')
    parser.add_argument('--add_module4', default='empty', help='Path to add_module4 file')
    parser.add_argument('--add_module5', default='empty', help='Path to add_module5 file')

    args = parser.parse_args()

    add_module_files = [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output, add_module_files)
