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

    # Append column names to genome_stats sheet
    column_names = ["sample", "number of scaffolds"]

    # Check if the columns exist in combined_annotations_df and append them if they do
    if "taxonomy" in combined_data.columns:
        column_names.append("taxonomy")
    if "Completeness" in combined_data.columns:
        column_names.append("completeness")
    if "Contamination" in combined_data.columns:
        column_names.append("contamination")

    column_names += ["tRNA count"]
    gs_sheet.append(column_names)

    # Populate genome_stats sheet with data from combined_annotations
    for sample in unique_samples:
        # Extract information for the current sample from combined_annotations
        sample_info = combined_data[combined_data['sample'] == sample]

        if not sample_info.empty:
            # Get the first row of sample_info (assuming one row per sample)
            sample_info = sample_info.iloc[0]

            # Extract completeness and contamination values
            completeness = sample_info.get('Completeness', None)
            contamination = sample_info.get('Contamination', None)
        else:
            # Handle the case where sample_info is empty (no data for the sample)
            completeness = None
            contamination = None

        # Append data to genome_stats sheet
        gs_data = [sample, None]

        # Check if the columns exist in combined_annotations_df and append them if they do
        if "taxonomy" in combined_data.columns:
            gs_data.append(sample_info.get('taxonomy', None))
        if "Completeness" in combined_data.columns:
            gs_data.append(completeness)
        if "Contamination" in combined_data.columns:
            gs_data.append(contamination)

        gs_data += [None]  # Placeholder for tRNA count, will be filled later
        gs_sheet.append(gs_data)

    # Sum tRNA counts and update the 'tRNA count' column
    trna_data = pd.read_csv(trna_file, sep='\t')

    for idx, sample in enumerate(unique_samples, start=2):
        # Extract relevant data for the current sample
        sample_trna_count = trna_data[sample].sum(numeric_only=True)

        # Update the 'tRNA count' column in the correct cell
        gs_sheet.cell(row=idx, column=len(column_names)).value = sample_trna_count

    # Create a dictionary to store data for each sheet
    sheet_data = {}

    for _, row in data.iterrows():
        # Split the "sheet" values by "; " and iterate over them
        for sheet_name in row['topic_ecosystem'].split('; '):
            sheet_name = sheet_name.replace(" ", "_")
            if sheet_name not in sheet_data:
                sheet_data[sheet_name] = {
                    'columns': [],  # Store column names for this sheet
                    'data': []  # Store data rows for this sheet
                }

            # Exclude the expected columns and move other columns
            row_data = [row[col] for col in row.index if col not in column_names]

            # Append the modified row to the corresponding sheet's data
            sheet_data[sheet_name]['data'].append(row_data)

            # Collect column names that are not in column_names
            new_columns = [col for col in row.index if col not in column_names]
            sheet_data[sheet_name]['columns'].extend(new_columns)

    # Inside the loop that creates sheets for topic_ecosystem values
    for sheet_name, sheet_info in sheet_data.items():
        # Create a worksheet for each sheet
        ws = wb.create_sheet(title=sheet_name)

        # Extract unique column names for this sheet based on the order of appearance
        sorted_column_names = list(sheet_info['columns'])

        # Append column names as the first row
        ws.append(sorted_column_names)

        # Append data rows to the worksheet
        for r_idx, row in enumerate(sheet_info['data'], 1):
            # Convert values in additional columns to strings
            row = [str(value) if col in sheet_info['columns'] else value for col, value in zip(sorted_column_names, row)]
            ws.append(row)

        # Create a table from the data for filtering
        tab = Table(displayName=f"{sheet_name}_Table", ref=ws.dimensions)
        style = TableStyleInfo(
            name="TableStyleMedium9", showFirstColumn=False,
            showLastColumn=False, showRowStripes=True, showColumnStripes=True
        )
        tab.tableStyleInfo = style
        ws.add_table(tab)

    # Before adding rRNA and tRNA sheets
    print("Adding rRNA sheet")
    # Add rRNA sheet
    rrna_data = pd.read_csv(rrna_file, sep='\t')
    rrna_sheet = wb.create_sheet(title="rRNA")

    # Append column names as the first row
    rrna_sheet.append(list(rrna_data.columns))

    # Append data rows to the worksheet
    for _, row in rrna_data.iterrows():
        rrna_sheet.append(list(row))

    # Before adding tRNA sheet
    print("Adding tRNA sheet")
    # Add tRNA sheet
    trna_data = pd.read_csv(trna_file, sep='\t')
    trna_sheet = wb.create_sheet(title="tRNA")

    # Append column names as the first row
    trna_sheet.append(list(trna_data.columns))

    # Append data rows to the worksheet
    for _, row in trna_data.iterrows():
        trna_sheet.append(list(row))

    # Before removing the default "Sheet" that was created
    print("Removing default 'Sheet'")
    # Remove the default "Sheet" that was created
    default_sheet = wb['Sheet']
    wb.remove(default_sheet)

    # Save the workbook as the output file
    print(f"Saving workbook to {output_file}")
    wb.save(output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate multi-sheet XLSX file from TSV files')
    parser.add_argument('--input_file', help='Input TSV file containing gene data')
    parser.add_argument('--rrna_file', help='rRNA TSV file')
    parser.add_argument('--trna_file', help='tRNA TSV file')
    parser.add_argument('--combined_annotations', help='Combined annotations TSV file')
    parser.add_argument('--combined_rrna', help='Combined rRNA TSV file')
    parser.add_argument('--output_file', help='Output XLSX file')
    args = parser.parse_args()

    generate_multi_sheet_xlsx(args.input_file, args.rrna_file, args.trna_file, args.combined_annotations, args.combined_rrna, args.output_file)
