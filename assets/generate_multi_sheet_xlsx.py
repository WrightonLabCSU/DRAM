import pandas as pd
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.worksheet.table import Table, TableStyleInfo

def generate_multi_sheet_xlsx(input_file, output_file):
    # Read the data from the input file using pandas with tab as the separator
    data = pd.read_csv(input_file, sep='\t')

    # Create a Workbook
    wb = Workbook()

    # Create a dictionary to store data for each sheet
    sheet_data = {}

    for _, row in data.iterrows():
        # Split the "sheet" values by "; " and iterate over them
        for sheet_name in row['sheet'].split('; '):
            sheet_name = sheet_name.replace(" ", "_")
            if sheet_name not in sheet_data:
                sheet_data[sheet_name] = []

            # Exclude the "sheet" column and move "gene_id" as the second column
            row_data = row.drop('sheet')
            gene_id = row['gene_id']
            row_data = pd.concat([pd.Series(gene_id), row_data[1:]], ignore_index=True)

            # Append the modified row to the corresponding sheet
            sheet_data[sheet_name].append(row_data)

    for sheet_name, sheet_rows in sheet_data.items():
        # Create a worksheet for each sheet
        ws = wb.create_sheet(title=sheet_name)

        # Extract column names from the original DataFrame
        column_names = data.columns.tolist()

        # Append column names as the first row
        ws.append(column_names)

        # Append data rows to the worksheet
        for r_idx, row in enumerate(dataframe_to_rows(pd.DataFrame(sheet_rows), index=False, header=False), 1):
            ws.append(row)

        # Create a table from the data for filtering
        tab = Table(displayName=f"{sheet_name}_Table", ref=ws.dimensions)
        style = TableStyleInfo(
            name="TableStyleMedium9", showFirstColumn=False,
            showLastColumn=False, showRowStripes=True, showColumnStripes=True
        )
        tab.tableStyleInfo = style
        ws.add_table(tab)

    # Remove the default "Sheet" that was created
    default_sheet = wb['Sheet']
    wb.remove(default_sheet)

    # Save the workbook as the output file
    wb.save(output_file)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate multi-sheet XLSX file')
    parser.add_argument('--input-file', required=True, help='Path to the input TSV file')
    parser.add_argument('--output-file', required=True, help='Path to the output XLSX file')

    args = parser.parse_args()
    generate_multi_sheet_xlsx(args.input_file, args.output_file)
