import pandas as pd
import sqlite3
from concurrent.futures import ThreadPoolExecutor
import argparse
import os

def fetch_descriptions(chunk, db_name, db_file):
    # Function to fetch descriptions based on IDs from the specified table
    table_name = f"{db_name}_description"
    ids_column = "id"  # Use the correct column name
    descriptions_column = f"{db_name}_description"
    
    # Establish connection to SQLite database
    conn = sqlite3.connect(db_file)
    
    ids = chunk[ids_column].unique()
    query = f"SELECT {ids_column}, {descriptions_column} FROM {table_name} WHERE {ids_column} IN ({','.join(['?'] * len(ids))})"
    
    cursor = conn.cursor()
    cursor.execute(query, ids)
    results = cursor.fetchall()
    
    descriptions_dict = {row[0]: row[1] for row in results}
    chunk[f"{db_name}_description"] = chunk[ids_column].map(descriptions_dict)
    
    # Close database connection
    conn.close()
    
    return chunk


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
    for i, chunk in enumerate(reader):
        if i >= 5:  # Print only the first 5 chunks
            break
        print("Chunk", i+1)
        print("Column Names:", chunk.columns)  # Print column names
        print(chunk.head())  # Print first few rows of the chunk

    # Process chunks concurrently
    with ThreadPoolExecutor() as executor:
        processed_chunks = executor.map(lambda chunk: fetch_descriptions(chunk, args.db_name, args.db_file), reader)

    # Concatenate processed chunks into a single DataFrame
    df = pd.concat(processed_chunks, ignore_index=True)

    # Write updated DataFrame to new CSV file
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
