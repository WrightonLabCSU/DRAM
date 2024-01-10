import pandas as pd
import argparse

def distill_summary(combined_annotations, genome_summary_form, target_id_counts, output, add_modules):
    # Read the input files into dataframes
    combined_data = pd.read_csv(combined_annotations, sep='\t')
    genome_summary = pd.read_csv(genome_summary_form, sep='\t')

    # Initialize a dictionary to store additional modules dataframes
    additional_modules = {}

    # Read and store data from add_moduleX files
    for i, add_module_file in enumerate(add_modules):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_modules[f'add_module{i + 1}'] = additional_module_data

    # Merge genome_summary and additional modules dataframes
    additional_modules_combined = pd.concat(additional_modules.values(), axis=1, keys=additional_modules.keys())

    # Identify the columns ending with "_id" for merging
    id_columns_combined = [col for col in combined_data.columns if col.endswith("_id")]

    # Exclude 'query_id' from the list of columns to be considered for merging
    id_columns_combined = [col for col in id_columns_combined if col != 'query_id']

    # Ensure 'gene_id' is included in the list of columns to be considered for merging
    id_columns_combined.append('gene_id')

    # Merge combined_data with genome_summary and additional modules based on columns ending with "_id"
    merged_data = pd.merge(
        combined_data,
        genome_summary,
        left_on=id_columns_combined,
        right_on=['gene_id'] * len(id_columns_combined),
        how='left'
    )

    # Create the distill_summary dataframe
    distill_summary = merged_data[['query_id', 'sample', 'gene_id'] + list(genome_summary.columns[1:])]

    # Write the distill_summary dataframe to the output file
    distill_summary.to_csv(output, sep='\t', index=False)

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
