import sqlite3
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Import annotations into an SQLite database.')
    parser.add_argument('--combined_annotations', type=str, help='Path to the combined annotations TSV file.')
    parser.add_argument('--db_name', type=str, help='Name of the SQLite database file.')
    parser.add_argument('--db_list', type=str, help='Space-separated list of database names used for annotation.')
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

def import_annotations(conn, file_path, db_list):
    df = pd.read_csv(file_path, sep='\t')
    db_names = db_list.split()  # Convert db_list string to a list of database names

    # Columns to extract for each database
    columns_to_extract = ['query_id', 'sample']
    for db_name in db_names:
        columns_to_extract += [col for col in df.columns if col.startswith(db_name) and ('_id' in col or '_EC' in col)]

    # Filter the DataFrame to only include the necessary columns
    df_filtered = df[columns_to_extract]

    # Melt the DataFrame to create separate rows for each ID/EC number
    df_melted = df_filtered.melt(id_vars=['query_id', 'sample'], var_name='db_column', value_name='gene_id').dropna()

    # Insert data into the database
    df_melted[['query_id', 'sample', 'gene_id']].to_sql('annotations', conn, if_exists='append', index=False, method='multi')

def main():
    args = parse_arguments()

    # Create the database and table
    conn = create_database(args.db_name)

    # Import annotations into the database
    import_annotations(conn, args.combined_annotations, args.db_list)

    conn.close()

if __name__ == '__main__':
    main()
