import pandas as pd
import sqlite3
from concurrent.futures import ThreadPoolExecutor
import argparse

def fetch_descriptions(chunk, db_name, db_file):
    # Function to fetch descriptions based on IDs from the specified table
    table_name = f"{db_name}_description"
    # Use "id" as the column name for fetching IDs from the hits CSV file
    ids_column = "id"
    descriptions_column = "description"
    
    # Establish connection to SQLite database
    conn = sqlite3.connect(db_file)
    
    # Adjust the column name to match the hits CSV file
    hits_ids_column = f"{db_name}_id"
    ids = chunk[hits_ids_column].unique()
    query = f"SELECT {ids_column}, {descriptions_column} FROM {table_name} WHERE {ids_column} IN ({','.join(['?'] * len(ids))})"
    
    cursor = conn.cursor()
    cursor.execute(query, ids)
    results = cursor.fetchall()
    
    descriptions_dict = {row[0]: row[1] for row in results}
    chunk[f"{db_name}_description"] = chunk[hits_ids_column].map(descriptions_dict)
    
    # Special processing for "kegg" database
    if db_name == "kegg":
        # Add additional output columns
        chunk["kegg_orthology"] = chunk[f"{db_name}_description"].apply(lambda x: extract_kegg_orthology(x))
        chunk["kegg_EC"] = chunk[f"{db_name}_description"].apply(lambda x: extract_kegg_EC(x))
    
    # Close database connection
    conn.close()
    
    return chunk


def extract_kegg_orthology(description):
    # Extract KO from the description
    if "(K" in description:
        ko_start = description.find("(K") + 1
        ko_end = description.find(")", ko_start)
        return description[ko_start:ko_end]
    else:
        return None


def extract_kegg_EC(description):
    # Extract EC numbers from the description
    ec_start = description.find("[EC:")
    if ec_start != -1:
        ec_end = description.find("]", ec_start)
        return description[ec_start + 5:ec_end]
    else:
        return None


def main():
    parser = argparse.ArgumentParser(description="Add descriptions from SQL database to hits file")
    parser.add_argument("--hits_csv", type=str, help="Path to the hits CSV file")
    parser.add_argument("--db_name", type=str, help="Name of the database table to fetch descriptions from")
    parser.add_argument("--output", type=str, help="Path to the output formatted CSV file")
    parser.add_argument("--db_file", type=str, help="Path to the SQLite database file")
    args = parser.parse_args()

    # Get the number of lines in the hits CSV file
    num_lines = sum(1 for line in open(args.hits_csv))

    # Determine chunk size dynamically based on number of lines
    chunksize = max(1, min(10000, num_lines // 10))  # Adjust as needed, minimum 1 chunk

    # Read CSV file in chunks
    reader = pd.read_csv(args.hits_csv, delimiter=',', chunksize=chunksize)

    # Process chunks
    with ThreadPoolExecutor() as executor:
        # Convert reader to list before passing to executor.map
        processed_chunks = executor.map(lambda chunk: fetch_descriptions(chunk, args.db_name, args.db_file), list(reader))

    # Concatenate processed chunks into a single DataFrame
    df = pd.concat(processed_chunks, ignore_index=True)

    # Write updated DataFrame to new CSV file
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
