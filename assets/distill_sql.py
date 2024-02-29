import sqlite3
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Import annotations into an SQLite database.')
    parser.add_argument('--combined_annotations', type=str, help='Path to the combined annotations TSV file.')
    parser.add_argument('--db_name', type=str, help='Name of the SQLite database file.')
    return parser.parse_args()

def create_database(db_name, include_extra_columns=False):
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    extra_columns_sql = ""
    if include_extra_columns:
        extra_columns_sql = ", taxonomy TEXT, Completeness REAL, Contamination REAL"
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS annotations (
            query_id TEXT,
            sample TEXT,
            gene_id TEXT{extra_columns_sql},
            UNIQUE(query_id, sample, gene_id)
        )
    ''')
    conn.commit()
    return conn

def identify_databases(df):
    db_names = set()
    for col in df.columns:
        if '_' in col:
            prefix = col.split('_')[0]
            db_names.add(prefix)
    non_db_prefixes = {'query', 'sample', 'start', 'end', 'strandedness', 'Completeness', 'Contamination', 'taxonomy'}
    return db_names - non_db_prefixes

def import_annotations(conn, file_path):
    df = pd.read_csv(file_path, sep='\t')
    db_names = identify_databases(df)

    extra_columns = [col for col in ['taxonomy', 'Completeness', 'Contamination'] if col in df.columns]

    data_to_insert = []
    for db_name in db_names:
        id_cols = [col for col in df.columns if col.startswith(db_name) and ('_id' in col or '_EC' in col)]
        for col in id_cols:
            for _, row in df[['query_id', 'sample', col] + extra_columns].dropna().iterrows():
                record = (row['query_id'], row['sample'], row[col]) + tuple(row[extra_column] for extra_column in extra_columns)
                data_to_insert.append(record)

    columns_sql = "query_id, sample, gene_id" + (", " + ", ".join(extra_columns) if extra_columns else "")
    placeholders_sql = ", ".join(["?"] * (3 + len(extra_columns)))
    cursor = conn.cursor()
    cursor.executemany(f'''
        INSERT OR IGNORE INTO annotations ({columns_sql})
        VALUES ({placeholders_sql})
    ''', data_to_insert)
    conn.commit()

def main():
    args = parse_arguments()
    df = pd.read_csv(args.combined_annotations, sep='\t')
    extra_columns = ['taxonomy', 'Completeness', 'Contamination']

    # Check if extra columns exist in the DataFrame
    include_extra_columns = any(col in df.columns for col in extra_columns)

    # Create the database and table
    conn = create_database(args.db_name, include_extra_columns=include_extra_columns)

    # Import annotations into the database
    import_annotations(conn, args.combined_annotations)

    conn.close()

if __name__ == '__main__':
    main()
