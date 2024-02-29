import pandas as pd
import sqlite3
import argparse
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate a multi-sheet XLSX document from distill sheets and a SQLite database.')
    parser.add_argument('--target_id_counts', type=str, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--db_name', type=str, help='Name of the SQLite database file.')
    parser.add_argument('--distill_sheets', nargs='+', help='List of paths to distill sheets.')
    parser.add_argument('--rrna_file', type=str, help='Path to the rrna_sheet.tsv file.', default=None)
    parser.add_argument('--trna_file', type=str, help='Path to the trna_sheet.tsv file.', default=None)
    parser.add_argument('--output_file', type=str, help='Path to the output XLSX file.')
    return parser.parse_args()

def compile_target_id_counts(target_id_counts):
    return pd.read_csv(target_id_counts, sep='\t')


def read_distill_sheets(distill_sheets):
    sheets_data = {}
    for sheet_path in distill_sheets:
        # Check for "NULL" content
        with open(sheet_path, 'r') as f:
            if f.readline().strip() == "NULL":
                print(f"Skipping {sheet_path} as it contains 'NULL'.")
                continue
        df = pd.read_csv(sheet_path, sep='\t')
        if 'topic_ecosystem' not in df.columns:
            print(f"Warning: 'topic_ecosystem' column not found in {sheet_path}. Skipping this file.")
            continue
        topic = df['topic_ecosystem'].unique().tolist()
        sheets_data.update({sheet_path: {'dataframe': df, 'topics': topic}})
    return sheets_data

def compile_genome_stats(db_name):
    conn = sqlite3.connect(db_name)
    query = "SELECT DISTINCT sample FROM annotations"
    df_samples = pd.read_sql_query(query, conn)
    # Assuming 'number_of_scaffolds' can be initially left blank or filled with placeholder values
    df_genome_stats = pd.DataFrame({
        "sample": df_samples['sample'],
        "number_of_scaffolds": ['' for _ in range(df_samples.shape[0])]  # Blank values
    })
    conn.close()
    return df_genome_stats


def add_sheet_from_dataframe(wb, df, sheet_name):
    ws = wb.create_sheet(title=sheet_name)
    for r in dataframe_to_rows(df, index=False, header=True):
        ws.append(r)

def main():
    args = parse_arguments()

    wb = Workbook()
    wb.remove(wb.active)  # Remove the default sheet

    # Compile and add genome stats sheet based on the database
    genome_stats_df = compile_genome_stats(args.db_name)
    add_sheet_from_dataframe(wb, genome_stats_df, "genome_stats")

    # Compile target_id_counts for later use
    target_id_counts_df = compile_target_id_counts(args.target_id_counts)

    # Read and process each distill sheet
    distill_data = read_distill_sheets(args.distill_sheets)
    for sheet_path, info in distill_data.items():
        for topic in info['topics']:
            # Your existing code to create topic-specific sheets and append data
            # Ensure to include the relevant target_id_counts for each gene_id in these sheets

    # Add tRNA and rRNA sheets if provided
    # Your existing code to handle tRNA and rRNA sheets

    wb.save(args.output_file)

if __name__ == '__main__':
    main()