import pandas as pd
import argparse

def distill_summary(combined_annotations, genome_summary_form, target_id_counts, output, add_modules):
    # Read input files
    combined_annotations_data = pd.read_csv(combined_annotations, sep='\t')
    genome_summary_form_data = pd.read_csv(genome_summary_form, sep='\t')

    # Select columns ending in '_id' for gene_id
    id_columns = [col for col in combined_annotations_data.columns if col.endswith('_id')]
    combined_annotations_data['gene_id'] = combined_annotations_data[id_columns].apply(lambda x: '; '.join(x.dropna()), axis=1)

    # Merge combined_annotations with genome_summary_form on gene_id
    merged_data = pd.merge(combined_annotations_data, genome_summary_form_data, left_on=['gene_id'], right_on=['gene_id'], how='left')

    # Create a dictionary to store additional modules data
    additional_modules = {}
    for i, add_module_file in enumerate(add_modules, start=1):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_modules[f'add_module{i}'] = additional_module_data

    # Select relevant columns from merged_data
    output_columns = ['query_id', 'sample'] + list(genome_summary_form_data.columns[1:])

    # Append additional modules data
    for key, additional_module_data in additional_modules.items():
        merged_data = pd.merge(merged_data, additional_module_data, on='gene_id', how='left', suffixes=('', f'_{key}'))

    # Write the final result to the output file
    merged_data[output_columns].to_csv(output, sep='\t', index=False)

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

    distill_summary(args.combined_annotations, args.genome_summary_form, args.target_id_counts, args.output, add_modules)
