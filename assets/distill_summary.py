import argparse
import pandas as pd

def distill_summary(combined_annotations_file, genome_summary_form_file, output_file, add_module_files=None):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_file, sep='\t')

    # Identify the column containing gene IDs dynamically
    gene_id_column = [col for col in combined_annotations.columns if col.endswith('_id') and col != 'query_id'][0]

    # Create a dictionary to store additional modules if provided
    additional_modules = {}
    if add_module_files:
        for i, add_module_file in enumerate(add_module_files):
            if add_module_file and add_module_file != 'empty':
                additional_module_data = pd.read_csv(add_module_file, sep='\t')
                additional_modules[f'add_module{i + 1}'] = additional_module_data

    # Merge genome_summary_form with combined_annotations based on gene_id
    merged_data = pd.merge(combined_annotations, genome_summary_form, left_on=gene_id_column, right_on='gene_id', how='inner')

    # Add additional modules if available
    for module_name, module_data in additional_modules.items():
        merged_data = pd.merge(merged_data, module_data, left_on=gene_id_column, right_on='gene_id', how='left')

    # Select columns for the distill_summary.tsv output
    output_columns = ['gene_id'] + list(genome_summary_form.columns[1:]) + list(additional_modules.keys())

    # Write the distilled summary to the output file
    merged_data[output_columns].to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate distill summary')
    parser.add_argument('--combined_annotations', required=True, help='Path to the combined annotations file')
    parser.add_argument('--genome_summary_form', required=True, help='Path to the genome summary file')
    parser.add_argument('--output', required=True, help='Path to the output file')
    parser.add_argument('--add_module1', required=False, help='Path to the additional module1 file')
    parser.add_argument('--add_module2', required=False, help='Path to the additional module2 file')
    parser.add_argument('--add_module3', required=False, help='Path to the additional module3 file')
    parser.add_argument('--add_module4', required=False, help='Path to the additional module4 file')
    parser.add_argument('--add_module5', required=False, help='Path to the additional module5 file')

    args = parser.parse_args()

    # Call the distill_summary function with the provided arguments
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output,
                    [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5])
