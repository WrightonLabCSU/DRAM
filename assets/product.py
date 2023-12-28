import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename='product.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def calculate_etc_modules(data):
    etc_modules_values = [None] * 20

    unique_module_names = data['module_name'].unique()
    for module_name in unique_module_names:
        module_data = data[data['module_name'] == module_name]
        total_annotations = float(module_data['ko'].count())

        for idx, sample in enumerate(module_data.columns[2:]):
            if module_data[sample].dtype != 'float64':
                continue
            sample_proportion = float(module_data[sample].sum()) / total_annotations
            etc_modules_values[idx] = sample_proportion

    return etc_modules_values

def create_product_tsv(target_counts_file, etc_file, module_step_file, function_file, output_file):
    logging.info("Loading target_id_counts data")
    target_counts_df = pd.read_csv(target_counts_file, sep='\t')
    
    logging.info("Loaded etc_module_database")
    etc_df = pd.read_csv(etc_file, sep='\t')
    
    logging.info("Loaded module_step_form")
    module_step_df = pd.read_csv(module_step_file, sep='\t')
    
    logging.info("Loaded function_heatmap_form")
    function_df = pd.read_csv(function_file, sep='\t')

    logging.info("Calculating ETC Modules values")
    etc_modules_values = calculate_etc_modules(module_step_df)

    # Create a new DataFrame for the final output
    output_data = module_step_df.copy()

    # Rename the 'ko' column to 'target_id' to match with target_counts_df
    output_data.rename(columns={'ko': 'target_id'}, inplace=True)

    # Set 'target_id' as the index for merging
    output_data.set_index('target_id', inplace=True)

    # Generate column mapping based on unique values in 'module_name'
    column_mapping = {col: f'Module_{idx}' for idx, col in enumerate(etc_df['module_name'].unique())}

    # Calculate proportions and add the results to the output_data DataFrame
    for col in output_data.columns:
        if output_data[col].dtype == 'float64':
            output_data[col] = output_data[col] / output_data[col].sum()

    # Reset the index
    output_data.reset_index(inplace=True)

    # Merge with target_counts_df using 'target_id'
    merged_data = target_counts_df.merge(output_data, on='target_id', how='left')

    # Create a new DataFrame for "ETC Modules"
    etc_modules_row = ["ETC Modules"] + [None] * (len(merged_data.columns) - 1)

    # Add the calculated etc_modules_values to the correct columns
    for idx, module_name in enumerate(etc_df['module_name'].unique()):
        if module_name in merged_data.columns:
            etc_modules_row[merged_data.columns.get_loc(module_name) + 1] = etc_modules_values[idx]

    # Add the "ETC Modules" row to the DataFrame
    merged_data.loc[len(merged_data)] = etc_modules_row

    # Save the output to a TSV file
    merged_data.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create product.tsv for heatmap')
    parser.add_argument('--input-target-counts', required=True, help='Path to target_id_counts file')
    parser.add_argument('--input-etc', required=True, help='Path to etc_module_database file')
    parser.add_argument('--input-module-step', required=True, help='Path to module_step_form file')
    parser.add_argument('--input-function-heatmap', required=True, help='Path to function_heatmap_form file')
    parser.add_argument('--output-file', required=True, help='Path to output product.tsv file')

    args = parser.parse_args()

    create_product_tsv(args.input_target_counts, args.input_etc, args.input_module_step, args.input_function_heatmap, args.output_file)
