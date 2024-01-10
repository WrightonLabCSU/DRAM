import pandas as pd
import argparse

def distill_summary(combined_annotations_file, genome_summary_form_file, target_id_counts_file, output_file, add_module_files):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_file, sep='\t')
    target_id_counts = pd.read_csv(target_id_counts_file, sep='\t')

    # Merge data from combined_annotations and target_id_counts based on sample
    merged_data = pd.merge(combined_annotations, target_id_counts, on='sample')

    # Create a dictionary to store additional module data
    additional_modules = {}
    for i, add_module_file in enumerate(add_module_files):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_modules[f'add_module{i+1}'] = additional_module_data

    # Process additional modules and merge with merged_data
    for module_name, module_data in additional_modules.items():
        merged_data = pd.merge(merged_data, module_data, on='gene_id', how='left')

        # Concatenate values if the column already exists in the output
        for column in module_data.columns[1:]:
            if column in merged_data.columns:
                merged_data[column] = merged_data[column].astype(str) + "; " + module_data[column]

    # Create the final output dataframe
    output_data = merged_data[['query_id', 'sample'] + list(combined_annotations.filter(regex='_id$').columns)]
    output_data = pd.merge(output_data, genome_summary_form, left_on='gene_id', right_on='gene_id', how='left')

    # Write the output to a TSV file
    output_data.to_csv(output_file, sep='\t', index=False)

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
