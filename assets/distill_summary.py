import argparse
import pandas as pd

def distill_summary(combined_annotations_path, genome_summary_form_path, output_path, add_module_files):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_path, sep='\t')

    # Additional columns from add_moduleX files
    for i, add_module_file in enumerate(add_module_files, start=1):
        if add_module_file != 'empty':
            add_module_data = pd.read_csv(add_module_file, sep='\t')

            # Rename columns to include add_moduleX prefix
            add_module_data = add_module_data.rename(lambda x: f'add_module{i}_{x}' if x != 'gene_id' else x, axis=1)

            # Merge add_moduleX data with genome_summary_form and concatenate values for duplicate columns
            genome_summary_form = pd.merge(genome_summary_form, add_module_data, on='gene_id', how='left')
            duplicate_cols = genome_summary_form.columns[genome_summary_form.columns.duplicated()]

            for col in duplicate_cols:
                if col.endswith('_sheet') or col.endswith('_module'):
                    # Concatenate values with a semicolon (;) for sheet and module columns
                    genome_summary_form[col[:-7]] = genome_summary_form[[col, f'{col}_add_module{i}']].apply(lambda x: '; '.join(x.dropna()), axis=1)
                else:
                    genome_summary_form[col] = genome_summary_form[[col, f'{col}_add_module{i}']].apply(lambda x: x.dropna().iloc[0] if x.dropna().any() else pd.NA, axis=1)

            # Drop the add_moduleX columns
            genome_summary_form = genome_summary_form.drop(columns=add_module_data.columns[1:])

    # Initialize the output DataFrame with gene_id, query_id, and sample columns
    distill_summary_df = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    # Additional columns from genome_summary_form
    additional_columns = list(genome_summary_form.columns)[1:]

    for col in additional_columns:
        distill_summary_df[col] = ''

    # Additional columns from combined_annotations (excluding "_id" columns and specific columns)
    combined_columns_to_add = [col for col in combined_annotations.columns
                               if not (col.endswith('_id') or col in ['query_id', 'banana_id', 'apple_id', 'pear_id', 'grape_id'])]

    for col in combined_columns_to_add:
        distill_summary_df[col] = ''

    # Iterate through gene_id values in genome_summary_form
    for _, row in genome_summary_form.iterrows():
        gene_id = row['gene_id']

        # Search for matches in combined_annotations columns ending in "_id"
        match_columns = [col for col in combined_annotations.columns if col.endswith('_id') and col != 'query_id']

        # Check if gene_id is present in combined_annotations
        if any(combined_annotations[column].eq(gene_id).any() for column in match_columns):
            for column in match_columns:
                # Check for matches
                match_rows = combined_annotations[combined_annotations[column] == gene_id]

                for _, match_row in match_rows.iterrows():
                    # Add matching information to the output DataFrame
                    distill_summary_df = distill_summary_df.append({
                        'gene_id': gene_id,
                        'query_id': match_row['query_id'],
                        'sample': match_row['sample'],
                    }, ignore_index=True)

                    # Add values from additional columns in genome_summary_form
                    for col in additional_columns:
                        distill_summary_df.at[distill_summary_df.index[-1], col] = row[col]

                    # Add values from selected columns in combined_annotations
                    for col in combined_columns_to_add:
                        distill_summary_df.at[distill_summary_df.index[-1], col] = match_row[col]

                    # Add values from add_moduleX columns, concatenate values with "; "
                    for add_module_col in add_module_data.columns[1:]:
                        new_col_name = add_module_col  # Keep the original column name without prefix

                        if new_col_name in distill_summary_df.columns:
                            existing_value = distill_summary_df.at[distill_summary_df.index[-1], new_col_name]
                            new_value = str(row[add_module_col])

                            # Check if the existing value is not NaN (numeric), then concatenate with "; "
                            if pd.notna(existing_value):
                                distill_summary_df.at[distill_summary_df.index[-1], new_col_name] = f'{existing_value}; {new_value}'
                            else:
                                distill_summary_df.at[distill_summary_df.index[-1], new_col_name] = new_value
                        else:
                            # Add new add_moduleX columns after existing columns
                            print(f"New Column Name: {new_col_name}")
                            print(f"Row Value: {row[add_module_col]}")
                            distill_summary_df[new_col_name] = str(row[add_module_col])


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

    # Call the distill_summary function with provided arguments
    distill_summary(
        args.combined_annotations,
        args.genome_summary_form,
        args.output,
        [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]
    )
