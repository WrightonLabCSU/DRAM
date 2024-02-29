import sqlite3
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Import annotations into an SQLite database.')
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
            UNIQUE(query_id, sample, gene_id)
        )
    ''')
    conn.commit()
    return conn

def identify_databases(df):
    # Extract database names from column headers
    db_names = set()
    for col in df.columns:
        if '_' in col:
            prefix = col.split('_')[0]
            # Assuming non-database columns don't follow the naming convention
            db_names.add(prefix)
    # Remove known non-database prefixes if any
    non_db_prefixes = {'query', 'sample', 'start', 'end', 'strandedness', 'Completeness', 'Contamination', 'taxonomy'}
    return db_names - non_db_prefixes

def import_annotations(conn, file_path):
    df = pd.read_csv(file_path, sep='\t')
    db_names = identify_databases(df)

    # Prepare data for insertion
    data_to_insert = []
    for db_name in db_names:
        id_cols = [col for col in df.columns if col.startswith(db_name) and ('_id' in col or '_EC' in col)]
        for col in id_cols:
            for _, row in df[['query_id', 'sample', col]].dropna().iterrows():
                data_to_insert.append((row['query_id'], row['sample'], row[col]))

    # Insert data into the database without duplicates
    cursor = conn.cursor()
    cursor.executemany('''
        INSERT OR IGNORE INTO annotations (query_id, sample, gene_id)
        VALUES (?, ?, ?)
    ''', data_to_insert)
    conn.commit()

def main():
    args = parse_arguments()

    # Create the database and table
    conn = create_database(args.db_name)

    # Import annotations into the database
    import_annotations(conn, args.combined_annotations)

    conn.close()

if __name__ == '__main__':
    main()
