import argparse
import pandas as pd
import os

def distill_summary(combined_annotations, genome_summary_form, output, add_modules):
    # Read input files
    combined_data = pd.read_csv(combined_annotations, sep='\t')
    genome_data = pd.read_csv(genome_summary_form, sep='\t')

    # Initialize an empty DataFrame for additional module data
    additional_modules = pd.DataFrame()

    # Read additional module files if provided
    for i, add_module_file in enumerate(add_modules, start=1):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_modules[f'add_module{i}'] = additional_module_data

    # Merge dataframes based on gene_id
    merged_data = pd.merge(combined_data, genome_data, left_on='gene_id', right_on='gene_id', how='inner')

    # Concatenate additional module data, if available
    if not additional_modules.empty:
        merged_data = pd.concat([merged_data, additional_modules], axis=1)

    # Check if output directory exists, create if not
    output_dir = os.path.dirname(output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write the result to the output file
    merged_data.to_csv(output, sep='\t', index=False)

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

    # Call the distill_summary function
    distill_summary(args.combined_annotations, args.genome_summary_form, args.output,
                    [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5])
