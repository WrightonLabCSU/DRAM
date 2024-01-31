import argparse
import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.worksheet.table import Table, TableStyleInfo

def generate_multi_sheet_xlsx(input_file, rrna_file, trna_file, output_file):
    # Read the data from the input file using pandas with tab as the separator
    data = pd.read_csv(input_file, sep='\t')

    # Create a Workbook
    wb = Workbook()

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
    parser.add_argument('--output-file', required=True, help='Path to the output XLSX file')

    args = parser.parse_args()
    generate_multi_sheet_xlsx(args.input_file, args.rrna_file, args.trna_file, args.output_file)
