import argparse
import pandas as pd

def distill_summary(combined_annotations_path, genome_summary_form_path, output_path, add_module_files):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_path, sep='\t')

    # Initialize the output DataFrame with gene_id, query_id, and sample columns
    distill_summary_df = pd.DataFrame(columns=['gene_id', 'query_id', 'sample'])

    # Additional columns from genome_summary_form
    additional_columns = ['gene_description', 'module', 'sheet', 'header', 'subheader']

    for col in additional_columns:
        distill_summary_df[col] = ''

    # Additional columns from combined_annotations (excluding "_id" columns and specific columns)
    combined_columns_to_add = [col for col in combined_annotations.columns
                               if not (col.endswith('_id') or col in ['query_id', 'banana_id', 'apple_id', 'pear_id', 'grape_id'])]

    for col in combined_columns_to_add:
        distill_summary_df[col] = ''

    # Iterate through gene_id values in genome_summary_form
    for gene_id in genome_summary_form['gene_id']:
        # Search for matches in combined_annotations columns ending in "_id"
        match_columns = [col for col in combined_annotations.columns if col.endswith('_id') and col != 'query_id']

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
                    # Concatenate values if there are multiple occurrences
                    values = genome_summary_form.loc[genome_summary_form['gene_id'] == gene_id, col].values
                    distill_summary_df.at[distill_summary_df.index[-1], col] = ';'.join(map(str, values))

                # Add values from selected columns in combined_annotations
                for col in combined_columns_to_add:
                    distill_summary_df.at[distill_summary_df.index[-1], col] = match_row[col]

    # Add values from additional module files
    for add_module_file in add_module_files:
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')

            for _, add_module_row in additional_module_data.iterrows():
                # Add matching information to the output DataFrame
                distill_summary_df = distill_summary_df.append({
                    'gene_id': add_module_row['gene_id'],
                    'query_id': add_module_row['query_id'],
                    'sample': add_module_row['sample'],
                }, ignore_index=True)

                # Add values from additional columns in additional module file
                for col in additional_columns:
                    # Concatenate values if there are multiple occurrences
                    values = additional_module_data.loc[
                        (additional_module_data['gene_id'] == add_module_row['gene_id']) &
                        (additional_module_data['query_id'] == add_module_row['query_id']) &
                        (additional_module_data['sample'] == add_module_row['sample']), col].values
                    distill_summary_df.at[distill_summary_df.index[-1], col] = ';'.join(map(str, values))

    # Write the output to a TSV file
    distill_summary_df.to_csv(output_path, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')
    parser.add_argument('--add_module1', default='empty', help='Path to the additional module file 1')
    parser.add_argument('--add_module2', default='empty', help='Path to the additional module file 2')
    parser.add_argument('--add_module3', default='empty', help='Path to the additional module file 3')
    parser.add_argument('--add_module4', default='empty', help='Path to the additional module file 4')
    parser.add_argument('--add_module5', default='empty', help='Path to the additional module file 5')

    args = parser.parse_args()

    add_module_files = [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]

    # Call the distill_summary function with provided arguments
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output, add_module_files)
