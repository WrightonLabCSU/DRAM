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
        if file_contains_data(sheet_path):
            df = pd.read_csv(sheet_path, sep='\t')
            topic = df['topic_ecosystem'].unique().tolist()
            sheets_data.update({sheet_path: {'dataframe': df, 'topics': topic}})
        else:
            print(f"Skipping {sheet_path} as it contains 'NULL'.")
    return sheets_data

def file_contains_data(file_path):
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            return first_line != "NULL"
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return False

def compile_genome_stats(db_name):
    conn = sqlite3.connect(db_name)
    df_genome_stats = pd.read_sql_query("SELECT sample, AVG(Completeness) AS Completeness, AVG(Contamination) AS Contamination, GROUP_CONCAT(DISTINCT taxonomy) AS taxonomy FROM annotations GROUP BY sample", conn)
    conn.close()
    return df_genome_stats

def query_annotations_for_gene_ids(db_name, gene_ids):
    conn = sqlite3.connect(db_name)
    placeholders = ', '.join('?' for _ in gene_ids)
    query = f"SELECT gene_id, sample FROM annotations WHERE gene_id IN ({placeholders})"
    df = pd.read_sql_query(query, conn, params=gene_ids)
    conn.close()
    return df

def add_sheet_from_dataframe(wb, df, sheet_name):
    ws = wb.create_sheet(title=sheet_name)
    for r in dataframe_to_rows(df, index=False, header=True):
        ws.append(r)

def main():
    args = parse_arguments()

    wb = Workbook()
    wb.remove(wb.active)  # Remove the default sheet

    genome_stats_df = compile_genome_stats(args.db_name)
    add_sheet_from_dataframe(wb, genome_stats_df, "genome_stats")

    target_id_counts_df = compile_target_id_counts(args.target_id_counts)
    distill_data = read_distill_sheets(args.distill_sheets)
    for sheet_path, info in distill_data.items():
        df_distill = info['dataframe']
        for topic in info['topics']:
            df_topic = df_distill[df_distill['topic_ecosystem'] == topic]
            gene_ids = df_topic['gene_id'].unique().tolist()
            df_annotations = query_annotations_for_gene_ids(args.db_name, gene_ids)
            # Merge directly with target_id_counts_df since the column names are consistent
            df_merged = pd.merge(df_topic, target_id_counts_df, on="gene_id", how="left")
            sheet_name = topic[:31]  # Excel sheet name character limit
            add_sheet_from_dataframe(wb, df_merged, sheet_name)

    # Add rrna and trna sheets if they contain data
    add_rrna_trna_sheets(wb, args.rrna_file, args.trna_file)

    wb.save(args.output_file)

def add_rrna_trna_sheets(wb, rrna_file, trna_file):
    for tsv_file, sheet_name in [(rrna_file, 'rRNA'), (trna_file, 'tRNA')]:
        if tsv_file and file_contains_data(tsv_file):
            df = pd.read_csv(tsv_file, sep='\t')
            add_sheet_from_dataframe(wb, df, sheet_name)

if __name__ == '__main__':
    main()
