import argparse
import pandas as pd
from openpyxl import Workbook
from openpyxl.worksheet.table import Table, TableStyleInfo

def generate_multi_sheet_xlsx(input_file, rrna_file, trna_file, combined_annotations, combined_rrna, output_file):
    # Read the data from the input file using pandas with tab as the separator
    data = pd.read_csv(input_file, sep='\t')

    # Read combined_annotations
    combined_data = pd.read_csv(combined_annotations, sep='\t')

    # Deduplicate and extract unique sample values
    unique_samples = combined_data['sample'].unique()

    # Create a Workbook
    wb = Workbook()

    # Create the genome_stats sheet
    gs_sheet = wb.create_sheet(title="genome_stats")

    # Dynamically get unique RNA types from combined_rrna
    rrna_data = pd.read_csv(combined_rrna, sep='\t')
    unique_rna_types = rrna_data['type'].unique()

    # Append column names to genome_stats sheet
    column_names = ["sample", "number of scaffolds", "taxonomy", "completeness", "contamination"] + list(unique_rna_types) + ["tRNA count"]
    gs_sheet.append(column_names)

    # Populate genome_stats sheet with data from combined_annotations
    for sample in unique_samples:
        sample_info = combined_data[combined_data['sample'] == sample].iloc[0]
        row_data = [
            sample,
            None,  # Placeholder for 'number of scaffolds'
            sample_info['taxonomy'],
            sample_info['Completeness'],
            sample_info['Contamination']
        ] + [None] * len(unique_rna_types) + [None]  # Placeholders for RNA counts and tRNA count
        gs_sheet.append(row_data)

    # Update RNA columns dynamically
    for rna_type in unique_rna_types:
        col_idx = column_names.index(rna_type) + 1
        for idx, sample in enumerate(unique_samples, start=2):
            sample_rrna_data = rrna_data[(rrna_data['sample'] == sample) & (rrna_data['type'] == rna_type)]
            if not sample_rrna_data.empty:
                values = [f"{row['query_id']}, ({row['begin']}, {row['end']})" for _, row in sample_rrna_data.iterrows()]
                joined_values = "; ".join(values)
                gs_sheet.cell(row=idx, column=col_idx).value = joined_values

    # Sum tRNA counts and update the 'tRNA count' column
    trna_data = pd.read_csv(trna_file, sep='\t')
    for idx, sample in enumerate(unique_samples, start=2):
        sample_trna_count = trna_data[trna_data['sample'] == sample]['count'].sum()
        gs_sheet.cell(row=idx, column=column_names.index("tRNA count") + 1).value = sample_trna_count

    # Create a dictionary to store data for each sheet
    sheet_data = {}

    # Fixed columns
    fixed_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']
    all_columns = fixed_columns.copy()
    if 'potential_amg' in data.columns and 'potential_amg' not in all_columns:
        all_columns.append('potential_amg')
    all_columns.extend(col for col in data.columns if col not in all_columns)

    for _, row in data.iterrows():
        for sheet_name in row['topic_ecosystem'].split('; '):
            sheet_name = sheet_name.replace(" ", "_")
            if sheet_name not in sheet_data:
                sheet_data[sheet_name] = []
            row_data = [row[col] for col in all_columns]
            sheet_data[sheet_name].append(row_data)

    for sheet_name, sheet_rows in sheet_data.items():
        ws = wb.create_sheet(title=sheet_name)
        ws.append(all_columns)
        for row in sheet_rows:
            ws.append(row)

    # Add rRNA sheet
    rrna_sheet = wb.create_sheet(title="rRNA")
    rrna_sheet.append(list(rrna_data.columns))
    for _, row in rrna_data.iterrows():
        rrna_sheet.append(list(row))

    # Remove the default "Sheet" that was created
    if 'Sheet' in wb.sheetnames:
        default_sheet = wb['Sheet']
        wb.remove(default_sheet)

    # Save the workbook as the output file
    wb.save(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate multi-sheet XLSX file')
    parser.add_argument('--input-file', required=True, help='Path to the input TSV file')
    parser.add_argument('--rrna-file', required=True, help='Path to the rRNA TSV file')
    parser.add_argument('--trna-file', required=True, help='Path to the tRNA TSV file')
    parser.add_argument('--combined-annotations', required=True, help='Path to the combined_annotations TSV file')
    parser.add_argument('--combined-rrna', required=True, help='Path to the combined_rrna TSV file')
    parser.add_argument('--output-file', required=True, help='Path to the output XLSX file')

    args = parser.parse_args()
    generate_multi_sheet_xlsx(args.input_file, args.rrna_file, args.trna_file, args.combined_annotations, args.combined_rrna, args.output_file)
