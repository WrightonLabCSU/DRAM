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

def file_contains_data(file_path):
    """Check if the file contains data other than 'NULL'."""
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            if first_line == "NULL":
                return False
        return True
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return False

def add_rrna_trna_sheets(wb, rrna_file, trna_file):
    """Read and add rRNA and tRNA sheets to the workbook if they contain data."""
    if rrna_file and file_contains_data(rrna_file):
        df_rrna = pd.read_csv(rrna_file, sep='\t')
        add_sheet_from_dataframe(wb, df_rrna, "rRNA")
    else:
        print(f"Skipping {rrna_file} as it contains 'NULL' or is not accessible.")

    if trna_file and file_contains_data(trna_file):
        df_trna = pd.read_csv(trna_file, sep='\t')
        add_sheet_from_dataframe(wb, df_trna, "tRNA")
    else:
        print(f"Skipping {trna_file} as it contains 'NULL' or is not accessible.")


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
    placeholders = ', '.join('?' for _ in gene_ids)
    # Alias gene_id as target_id in the query
    query = f"SELECT gene_id AS target_id, sample FROM annotations WHERE gene_id IN ({placeholders})"
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

    # Add genome_stats sheet
    genome_stats_df = compile_genome_stats(args.db_name)
    add_sheet_from_dataframe(wb, genome_stats_df, "genome_stats")

    # Process and add distill sheets
    target_id_counts_df = compile_target_id_counts(args.target_id_counts)
    distill_data = read_distill_sheets(args.distill_sheets)
    for sheet_path, info in distill_data.items():
        df_distill = info['dataframe']
        for topic in info['topics']:
            df_topic = df_distill[df_distill['topic_ecosystem'] == topic]
            gene_ids = df_topic['gene_id'].unique().tolist()
            df_annotations = query_annotations_for_gene_ids(args.db_name, gene_ids)
            df_annotations.rename(columns={'gene_id': 'target_id'}, inplace=True)
            df_merged = pd.merge(df_topic, df_annotations, on="target_id", how="left")
            df_merged = pd.merge(df_merged, target_id_counts_df, on="target_id", how="left")
            df_final = df_merged.drop(columns=['query_id', 'sample', 'taxonomy', 'Completeness', 'Contamination'], errors='ignore')
            sheet_name = topic[:31]  # Excel sheet name character limit
            add_sheet_from_dataframe(wb, df_final, sheet_name)

    # Add rrna and trna sheets if not NULL
    for tsv_file, sheet_name in [(args.rrna_file, 'rRNA'), (args.trna_file, 'tRNA')]:
        if tsv_file:  # Check if file path is provided
            with open(tsv_file, 'r') as f:
                if f.readline().strip() != "NULL":  # File has content other than "NULL"
                    f.seek(0)  # Reset file pointer to the beginning of the file
                    df = pd.read_csv(f, sep='\t')
                    add_sheet_from_dataframe(wb, df, sheet_name)
                else:
                    print(f"Skipping {tsv_file} as it contains 'NULL'.")

    # Save the workbook
    wb.save(args.output_file)

if __name__ == '__main__':
    main()


if __name__ == '__main__':
    main()

