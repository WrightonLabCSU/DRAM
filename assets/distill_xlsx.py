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
    cursor = conn.cursor()
    
    # Check for the presence of optional columns in the database
    cursor.execute("PRAGMA table_info(annotations)")
    columns_info = cursor.fetchall()
    column_names = [info[1] for info in columns_info]
    
    # Base query to select distinct samples
    query = "SELECT DISTINCT sample FROM annotations"
    df_samples = pd.read_sql_query(query, conn)
    
    # Initialize the genome_stats DataFrame
    df_genome_stats = pd.DataFrame({
        "sample": df_samples['sample']
    })
    
    # If optional columns are present, aggregate and add them to df_genome_stats
    optional_columns = ['taxonomy', 'Completeness', 'Contamination']
    for col in optional_columns:
        if col in column_names:
            query = f"SELECT sample, AVG({col}) AS {col} FROM annotations GROUP BY sample"
            df_col_stats = pd.read_sql_query(query, conn)
            df_genome_stats = pd.merge(df_genome_stats, df_col_stats, on="sample", how="left")
    
    conn.close()
    return df_genome_stats


def query_annotations_for_gene_ids(db_name, gene_ids):
    conn = sqlite3.connect(db_name)
    placeholders = ', '.join('?' for gene_id in gene_ids)
    query = f"SELECT * FROM annotations WHERE gene_id IN ({placeholders})"
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
            # Rename 'gene_id' to 'target_id' in df_annotations for consistency
            df_annotations.rename(columns={'gene_id': 'target_id'}, inplace=True)

            df_topic_renamed = df_topic.rename(columns={'gene_id': 'target_id'})
            # Now both DataFrames have the 'target_id' column for a consistent merge
            df_merged = pd.merge(df_topic_renamed, df_annotations, on="target_id", how="left")
            df_merged = pd.merge(df_merged, target_id_counts_df, on="target_id", how="left")

            sheet_name = topic[:31]  # Excel sheet name character limit
            add_sheet_from_dataframe(wb, df_merged, sheet_name)

    wb.save(args.output_file)

if __name__ == '__main__':
    main()

