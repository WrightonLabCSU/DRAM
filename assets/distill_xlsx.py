import pandas as pd
import sqlite3
import argparse
from openpyxl import Workbook

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate a multi-sheet XLSX document from distill sheets and a SQLite database.')
    parser.add_argument('--target_id_counts', type=str, help='Path to the target_id_counts.tsv file.')
    parser.add_argument('--db_name', type=str, help='Name of the SQLite database file.')
    parser.add_argument('--distill_sheets', nargs='+', help='List of paths to distill sheets.')
    parser.add_argument('--rrna_file', type=str, help='Path to the rrna_sheet.tsv file.', default=None)
    parser.add_argument('--trna_file', type=str, help='Path to the trna_sheet.tsv file.', default=None)
    parser.add_argument('--output_file', type=str, help='Path to the output XLSX file.')
    return parser.parse_args()

def read_distill_sheets(distill_sheets):
    sheets_data = {}
    for sheet_path in distill_sheets:
        # Open the file and check if its contents are just "NULL"
        with open(sheet_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line == "NULL":
                print(f"Skipping {sheet_path} as it contains 'NULL'.")
                continue  # Skip this file and move to the next

        # If the file is not skipped, proceed to read it into a DataFrame
        df = pd.read_csv(sheet_path, sep='\t')
        if 'topic_ecosystem' not in df.columns:
            print(f"Warning: 'topic_ecosystem' column not found in {sheet_path}. Skipping this file.")
            continue  # Skip this file and move to the next
        topic = df['topic_ecosystem'].unique().tolist()
        sheets_data.update({sheet_path: {'dataframe': df, 'topics': topic}})
    return sheets_data


def query_database_for_gene_ids(db_name, gene_ids):
    conn = sqlite3.connect(db_name)
    placeholder= '?' # For SQLite. Adjust for other databases.
    placeholders= ', '.join(placeholder for unused in gene_ids)
    query = f'SELECT query_id, sample, gene_id FROM annotations WHERE gene_id IN ({placeholders})'
    df = pd.read_sql_query(query, conn, params=gene_ids)
    conn.close()
    return df

def compile_genome_stats(target_id_counts):
    return pd.read_csv(target_id_counts, sep='\t')

def create_topic_sheets(wb, distill_data, db_name):
    for sheet_path, info in distill_data.items():
        for topic in info['topics']:
            # Query for each gene_id within this topic
            gene_ids = info['dataframe']['gene_id'].tolist()
            matched_annotations = query_database_for_gene_ids(db_name, gene_ids)
            # Here, you would merge/join `matched_annotations` with `info['dataframe']` based on `gene_id`
            # and filter by `topic`, then add to the workbook

def add_trna_rrna_sheets(wb, rrna_file, trna_file):
    if rrna_file:
        df_rrna = pd.read_csv(rrna_file, sep='\t')
        # Add rrna sheet to workbook
    if trna_file:
        df_trna = pd.read_csv(trna_file, sep='\t')
        # Add trna sheet to workbook

def main():
    args = parse_arguments()

    wb = Workbook()
    distill_data = read_distill_sheets(args.distill_sheets)

    # Compile and add genome stats sheet
    genome_stats_df = compile_genome_stats(args.target_id_counts)
    # Add genome stats sheet to workbook

    create_topic_sheets(wb, distill_data, args.db_name)

    add_trna_rrna_sheets(wb, args.rrna_file, args.trna_file)

    wb.save(args.output_file)

if __name__ == '__main__':
    main()
