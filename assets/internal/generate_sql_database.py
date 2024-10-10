import os
import sqlite3
import argparse

def insert_data(conn, table_name, data):
    placeholders = ', '.join(['?'] * len(data[0]))
    query = f"INSERT OR REPLACE INTO {table_name} VALUES ({placeholders})"
    conn.executemany(query, data)
    conn.commit()

def process_dbcan(db_dir):
    description_file = os.path.join(db_dir, 'dbcan.fam-activities.tsv')
    ec_file = os.path.join(db_dir, 'dbcan.fam.subfam.ec.tsv')

    descriptions = {}
    ecs = {}
    skipped_lines = []

    # Process descriptions
    with open(description_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                descriptions[parts[0]] = ' '.join(parts[1:])
            elif len(parts) == 1:
                descriptions[parts[0]] = "No description available"
            else:
                skipped_lines.append(f"Skipped line in description file: {line.strip()} (expected at least 2 columns, found {len(parts)})")

    # Process EC numbers
    with open(ec_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 2:
                ecs[parts[0]] = ecs.get(parts[0], set())
                ecs[parts[0]].add(parts[2])

    data = []
    for entry in descriptions:
        ec = ','.join(ecs.get(entry, []))
        data.append((entry, descriptions[entry], ec))
    
    return data, skipped_lines

def main():
    parser = argparse.ArgumentParser(description="Generate descriptions database for DRAM2.")
    parser.add_argument('--db_dir', required=True, help="Directory containing the database subdirectories.")
    parser.add_argument('--output_db', required=True, help="Path to the output SQLite database.")
    parser.add_argument('--log', required=True, help="Path to the log file.")

    args = parser.parse_args()
    
    log_entries = []
    db_dir = args.db_dir
    output_db = args.output_db

    conn = sqlite3.connect(output_db)
    log_entries.append(f"Opened database {output_db}")

    dbcan_dir = os.path.join(db_dir, 'dbcan')
    if os.path.exists(dbcan_dir):
        conn.execute("""
            CREATE TABLE IF NOT EXISTS dbcan_description (
                id VARCHAR(30) NOT NULL, 
                description VARCHAR(1000), 
                ec VARCHAR(1000), 
                PRIMARY KEY (id)
            );
        """)
        log_entries.append("Processing dbcan_description from " + dbcan_dir)
        data, skipped_lines = process_dbcan(dbcan_dir)
        insert_data(conn, 'dbcan_description', data)
        log_entries.append(f"Inserted {len(data)} records into dbcan_description")
        log_entries.extend(skipped_lines)
    
    with open(args.log, 'w') as log_file:
        for entry in log_entries:
            log_file.write(entry + '\n')
    
    conn.close()
    log_entries.append("Closed database connection")

if __name__ == "__main__":
    main()