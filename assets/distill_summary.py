import pandas as pd
import argparse

def distill_summary(combined_annotations, genome_summary_form, target_id_counts, output, add_modules):
    # Read input files
    combined_annotations_data = pd.read_csv(combined_annotations, sep='\t')
    genome_summary_form_data = pd.read_csv(genome_summary_form, sep='\t')

    # Create a dictionary to store additional modules data
    additional_modules = {}
    for i, add_module_file in enumerate(add_modules, start=1):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_modules[f'add_module{i}'] = additional_module_data

    # Select columns ending in '_id' for gene_id
    id_columns = [col for col in combined_annotations_data.columns if col.endswith('_id')]
    combined_annotations_data['gene_id'] = combined_annotations_data[id_columns].apply(lambda x: '; '.join(x.dropna()), axis=1)

    # Merge combined_annotations with genome_summary_form on gene_id
    merged_data = pd.merge(combined_annotations_data, genome_summary_form_data, on='gene_id', how='inner')

    # Select relevant columns from merged_data
    output_columns = ['query_id', 'sample', 'gene_id'] + list(genome_summary_form_data.columns[1:])

    # Create a DataFrame to store the final result
    final_data = pd.DataFrame(columns=output_columns)

    # Iterate over rows in merged_data to handle multiple gene_id values
    for _, row in merged_data.iterrows():
        gene_ids = row['gene_id'].split('; ')
        for gene_id in gene_ids:
            final_data = final_data.append(row[output_columns[:-1]] + pd.Series([gene_id]), ignore_index=True)

    # Merge with target_id_counts on 'sample' column
    final_data = pd.merge(final_data, target_id_counts, left_on='sample', right_on='sample', how='left').drop(columns='gene_id_y').rename(columns={'gene_id_x': 'gene_id'})

    # Append additional modules data
    for key, additional_module_data in additional_modules.items():
        final_data = pd.merge(final_data, additional_module_data, on='gene_id', how='left', suffixes=('', f'_{key}'))

    # Write the final result to the output file
    final_data.to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--target_id_counts', required=True, help='Path to the target_id_counts file')
    parser.add_argument('--output', required=True, help='Path to the output file')
    parser.add_argument('--add_module1', required=False, help='Path to the additional module1 file')
    parser.add_argument('--add_module2', required=False, help='Path to the additional module2 file')
    parser.add_argument('--add_module3', required=False, help='Path to the additional module3 file')
    parser.add_argument('--add_module4', required=False, help='Path to the additional module4 file')
    parser.add_argument('--add_module5', required=False, help='Path to the additional module5 file')

    args = parser.parse_args()
    add_modules = [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]

    target_id_counts_data = pd.read_csv(args.target_id_counts, sep='\t')

    distill_summary(args.combined_annotations, args.genome_summary_form, target_id_counts_data, args.output, add_modules)
