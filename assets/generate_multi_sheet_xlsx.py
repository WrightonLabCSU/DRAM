import argparse
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
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

    # Append column names to genome_stats sheet
    gs_sheet.append(["sample", "number of scaffolds", "taxonomy", "completeness", "contamination", "5S rRNA", "16S rRNA", "23S rRNA", "tRNA count"])

    # Populate genome_stats sheet with data from combined_annotations
    for sample in unique_samples:
        # Extract information for the current sample from combined_annotations
        sample_info = combined_data[combined_data['sample'] == sample].iloc[0]  # Assuming one row per sample

        # Append data to genome_stats sheet
        gs_sheet.append([sample, None, sample_info.get('taxonomy', None), sample_info.get('Completeness', None), sample_info.get('Contamination', None), None, None, None, None])

    # Read combined_rrna
    rrna_data = pd.read_csv(combined_rrna, sep='\t')

    for rna_type in ["5S rRNA", "16S rRNA", "23S rRNA"]:
        # Filter rrna_data based on rna_type
        filtered_rrna_data = rrna_data[rrna_data['type'] == rna_type]

        print(f"RNA Type: {rna_type}")
        print(filtered_rrna_data)

        # Iterate over unique samples
        for sample in unique_samples:
            # Extract relevant data for the current sample and rna_type
            sample_rrna_data = filtered_rrna_data[filtered_rrna_data['sample'] == sample]

            print(f"Sample: {sample}")
            print(sample_rrna_data)

            # Concatenate the values and format them as needed
            values = [
                f"{row['query_id']}, ({row['begin']}, {row['end']})"
                for _, row in sample_rrna_data.iterrows()
            ]

            # Join multiple values with "; "
            joined_values = "; ".join(values)

            # Find the corresponding column index in the genome_stats sheet and update the value
            for col_idx, col in enumerate(gs_sheet[1], start=1):
                if col.value == rna_type:
                    for row in gs_sheet.iter_rows(min_row=2, max_row=gs_sheet.max_row, min_col=col_idx, max_col=col_idx):
                        if row[0].value == sample:
                            row_idx = row[0].row
                            gs_sheet.cell(row=row_idx, column=col_idx).value = joined_values


    # Create a dictionary to store data for each sheet
    sheet_data = {}

    # Fixed columns
    fixed_columns = ['gene_id', 'gene_description', 'pathway', 'topic_ecosystem', 'category', 'subcategory']

    for _, row in data.iterrows():
        # Split the "sheet" values by "; " and iterate over them
        for sheet_name in row['topic_ecosystem'].split('; '):  # Assuming 'topic_ecosystem' corresponds to 'sheet'
            sheet_name = sheet_name.replace(" ", "_")
            if sheet_name not in sheet_data:
                sheet_data[sheet_name] = []

            # Exclude the "sheet" column and move "gene_id" as the second column
            row_data = [row[col] for col in fixed_columns]

            # Include the 'potential_amg' column if it exists
            if 'potential_amg' in data.columns:
                # Convert 'potential_amg' values to "TRUE" or "FALSE"
                row_data += ['TRUE' if row['potential_amg'] == 'TRUE' else 'FALSE']

            # Append the rest of the columns without 'potential_amg'
            row_data += [row[col] for col in data.columns if col not in fixed_columns and col != 'potential_amg']

            # Append the modified row to the corresponding sheet
            sheet_data[sheet_name].append(row_data)

    for sheet_name, sheet_rows in sheet_data.items():
        # Create a worksheet for each sheet
        ws = wb.create_sheet(title=sheet_name)

        # Extract column names from the original DataFrame
        column_names = fixed_columns

        # Include 'potential_amg' column if it exists
        if 'potential_amg' in data.columns:
            column_names += ['potential_amg']

        # Append the rest of the columns without 'potential_amg'
        column_names += [col for col in data.columns if col not in fixed_columns and col != 'potential_amg']

        # Append column names as the first row
        ws.append(column_names)

        # Append data rows to the worksheet
        for r_idx, row in enumerate(sheet_rows, 1):
            ws.append(row)

        # Create a table from the data for filtering
        tab = Table(displayName=f"{sheet_name}_Table", ref=ws.dimensions)
        style = TableStyleInfo(
            name="TableStyleMedium9", showFirstColumn=False,
            showLastColumn=False, showRowStripes=True, showColumnStripes=True
        )
        tab.tableStyleInfo = style
        ws.add_table(tab)

    # Add rRNA and tRNA sheets
    rrna_data = pd.read_csv(rrna_file, sep='\t')
    trna_data = pd.read_csv(trna_file, sep='\t')

    rrna_sheet = wb.create_sheet(title="rRNA")
    trna_sheet = wb.create_sheet(title="tRNA")

    for sheet, data in zip([rrna_sheet, trna_sheet], [rrna_data, trna_data]):
        # Append column names as the first row
        sheet.append(list(data.columns))

        # Append data rows to the worksheet
        for _, row in data.iterrows():
            sheet.append(list(row))

    # Remove the default "Sheet" that was created
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
