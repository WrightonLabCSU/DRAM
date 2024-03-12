import pandas as pd
import sqlite3
import re
from concurrent.futures import ThreadPoolExecutor
import argparse

def fetch_descriptions(chunk, db_name, db_file):
    # Function to fetch descriptions based on IDs from the specified table
    table_name = f"{db_name}_description"
    ids_column = "id"  # Column name for fetching IDs
    descriptions_column = "description"
    
    # Establish connection to SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Adjust the column name to match the hits CSV file and process dbcan_id values if needed
    hits_ids_column = f"{db_name}_id"
    if db_name == "dbcan":
        chunk[hits_ids_column] = chunk[hits_ids_column].str.replace(".hmm", "")
    ids = chunk[hits_ids_column].unique()
    
    # Construct and execute the query
    if db_name == "dbcan":
        query = f"SELECT {ids_column}, {descriptions_column}, ec FROM {table_name} WHERE {ids_column} IN ({','.join(['?'] * len(ids))})"
    else:
        query = f"SELECT {ids_column}, {descriptions_column} FROM {table_name} WHERE {ids_column} IN ({','.join(['?'] * len(ids))})"
    
    cursor.execute(query, ids)
    results = cursor.fetchall()
    
    # Map fetched descriptions to the chunk
    if db_name == "dbcan":
        descriptions_dict = {row[0]: (row[1], row[2]) for row in results}
    else:
        descriptions_dict = {row[0]: row[1] for row in results}
    
    chunk[f"{db_name}_description"] = chunk[hits_ids_column].map(lambda x: descriptions_dict.get(x, ("", ""))[0])
    
    # Directly use the 'ec' column for dbcan or format and extract EC numbers from description for kegg
    if db_name == "dbcan":
        # Format the EC numbers with "EC:" prefix and semicolon-separated for dbcan
        chunk["dbcan_EC"] = chunk[hits_ids_column].map(lambda x: '; '.join([f"EC:{ec}" for ec in descriptions_dict.get(x, ("", ""))[1].split(',') if ec]) if descriptions_dict.get(x, ("", ""))[1] else "")
    elif db_name == "kegg":
        # Use the original description text to extract and format EC numbers for kegg
        chunk["kegg_EC"] = chunk[f"{db_name}_description"].apply(format_and_extract_EC_numbers)
    
    # Extract KEGG orthology if needed (for KEGG database)
    if db_name == "kegg":
        chunk["kegg_orthology"] = chunk[f"{db_name}_description"].apply(extract_kegg_orthology)
    
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


def format_and_extract_EC_numbers(description):
    # Find the start of the EC number section in the description
    ec_start = description.find("[EC:")
    if ec_start != -1:
        ec_end = description.find("]", ec_start)
        ec_text = description[ec_start + 4:ec_end]  # Skip the "[EC:" part
        # Extract EC numbers using regex and format them with "EC:" prefix
        ec_numbers = re.findall(r'\b\d+\.\d+\.\d+\.\d+\b', ec_text)
        # Format each EC number with "EC:" prefix and join with "; "
        formatted_ec_numbers = '; '.join([f"EC:{ec}" for ec in ec_numbers])
        return formatted_ec_numbers


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
