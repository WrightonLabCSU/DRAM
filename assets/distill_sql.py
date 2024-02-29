import sqlite3
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Import annotations into an SQLite database, dynamically identifying database columns.')
    parser.add_argument('--combined_annotations', type=str, help='Path to the combined annotations TSV file.')
    parser.add_argument('--db_name', type=str, help='Name of the SQLite database file.')
    return parser.parse_args()

def create_database(db_name):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS annotations (
            query_id TEXT,
            sample TEXT,
            gene_id TEXT,
            PRIMARY KEY (query_id, sample, gene_id)
        )
    ''')
    conn.commit()
    return conn

def identify_databases(column_names):
    # Define a list of known non-database column names
    non_db_columns = ['query_id', 'sample', 'start_position', 'end_position', 'strandedness', 'Completeness', 'Contamination', 'taxonomy']
    # Define common suffixes that might indicate a column is related to database annotations
    db_suffixes = ['_id', '_bitScore', '_rank', '_score_type', '_definition', '_EC', '_family', '_subfam_GenBank', '_subfam_EC']

    db_prefixes = set()
    for col in column_names:
        if any(col.endswith(suffix) for suffix in db_suffixes) and col not in non_db_columns:
            # Extract the database prefix by removing the suffix
            db_prefix = col.split('_')[0]
            db_prefixes.add(db_prefix)
    
    return list(db_prefixes)

def import_annotations(conn, file_path):
    df = pd.read_csv(file_path, sep='\t')
    db_names = identify_databases(df.columns)

    # Columns to extract for each database
    columns_to_extract = ['query_id', 'sample']
    for db_name in db_names:
        columns_to_extract += [col for col in df.columns if col.startswith(db_name)]

    # Filter the DataFrame to only include the necessary columns
    df_filtered = df[columns_to_extract]

    # For each database, extract ID and possibly EC numbers, and insert into the database
    for db_name in db_names:
        db_columns = [col for col in df_filtered.columns if col.startswith(db_name) and ('_id' in col or '_EC' in col)]
        for col in db_columns:
            # Create a temporary DataFrame for each ID/EC column to insert
            temp_df = df_filtered[['query_id', 'sample', col]].rename(columns={col: 'gene_id'}).dropna()

            # Drop duplicate rows based on 'query_id', 'sample', and 'gene_id' to avoid UNIQUE constraint violations
            temp_df = temp_df.drop_duplicates(subset=['query_id', 'sample', 'gene_id'])

            temp_df.to_sql('annotations', conn, if_exists='append', index=False, method='multi')


def main():
    args = parse_arguments()

    # Create the database and table
    conn = create_database(args.db_name)

    # Import annotations into the database, dynamically identifying database columns
    import_annotations(conn, args.combined_annotations)

    conn.close()

if __name__ == '__main__':
    main()
