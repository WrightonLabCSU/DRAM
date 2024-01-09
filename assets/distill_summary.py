import argparse
import pandas as pd
import logging

# Configure the logger
logging.basicConfig(filename="distill_summary.log", level=logging.INFO, format='%(levelname)s: %(message)s')

def expand_target_ids(data, target_id_col):
    rows = []
    for _, row in data.iterrows():
        target_ids = row[target_id_col].split('; ')
        for target_id in target_ids:
            new_row = row.copy()
            new_row[target_id_col] = target_id
            rows.append(new_row)
    return pd.DataFrame(rows)

def combine_dataframes(primary_df, additional_df):
    # Merge primary_df and additional_df using 'gene_id' column
    merged_df = primary_df.merge(additional_df, on='gene_id', how='outer')

    # Handle column collisions by concatenating values with "; " when both values are not empty
    for col in primary_df.columns:
        if col in additional_df.columns and col != 'gene_id':
            col_x = col + '_x'
            col_y = col + '_y'
            merged_df[col] = merged_df.apply(lambda row: f"{row[col_x]}; {row[col_y]}" if (not pd.isna(row[col_x]) and not pd.isna(row[col_y])) else row[col_x] if not pd.isna(row[col_x]) else row[col_y] if not pd.isna(row[col_y]) else '', axis=1)
            merged_df = merged_df.drop(columns=[col_x, col_y])

    return merged_df

def distill_summary(combined_annotations, genome_summary_form, target_id_counts, output_file, add_modules, target_id_col):
    # Log the input arguments
    logging.info(f"Combined Annotations file: {combined_annotations}")
    logging.info(f"Genome Summary Form file: {genome_summary_form}")
    logging.info(f"Target ID Counts file: {target_id_counts}")
    logging.info(f"Output file: {output_file}")
    logging.info(f"Add Modules: {add_modules}")
    logging.info(f"Target ID Column: {target_id_col}")

    # Read the data from the input files using pandas
    combined_annotations_data = pd.read_csv(combined_annotations, sep='\t')
    genome_summary_data = pd.read_csv(genome_summary_form, sep='\t')
    target_id_counts_data = pd.read_csv(target_id_counts, sep='\t')

    # Create an empty dictionary to store additional modules
    additional_modules = {}

    for i, add_module_file in enumerate(add_modules, 1):
        if add_module_file and add_module_file != 'empty':
            additional_module_data = pd.read_csv(add_module_file, sep='\t')
            additional_modules[f'add_module{i}'] = additional_module_data

    if additional_modules:
        # Start with the genome_summary_data
        combined_data = genome_summary_data

        # Merge with additional modules
        for module_name, additional_module_data in additional_modules.items():
            combined_data = combine_dataframes(combined_data, additional_module_data)

        # Merge the combined data with target_id_counts_data
        expanded_annotations = expand_target_ids(combined_annotations_data, target_id_col)
        merged_data = expanded_annotations.merge(combined_data, left_on=target_id_col, right_on='gene_id', how='inner')

        # Merge with target_id_counts, handle duplicate column names, and fill missing values with zeros
        combined_data = merged_data.merge(target_id_counts_data, left_on=target_id_col, right_on='target_id', how='left')
        combined_data = combined_data.fillna(0)

        # Remove the target_id column
        combined_data = combined_data.drop(columns=[target_id_col])

        # Log the number of rows in the output data
        logging.info(f"Number of rows in the output data: {len(combined_data)}")

        # Save the resulting dataframe to the output file in TSV format
        combined_data.to_csv(output_file, sep='\t', index=False)
        logging.info("Distill summary completed.")

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
    parser.add_argument('--target_id_col', required=True, help='Column name containing target ids in combined_annotations file')

    args = parser.parse_args()
    add_modules = [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]
    
    distill_summary(args.combined_annotations, args.genome_summary_form, args.target_id_counts, args.output, add_modules, args.target_id_col)
