import pandas as pd
import sqlite3
import re
from concurrent.futures import ThreadPoolExecutor
import argparse

def fetch_descriptions(chunk, db_name, db_file):
    # Function to fetch descriptions based on IDs from the specified table
    table_name = f"{db_name}_description"
    ids_column = "id"
    descriptions_column = "description"
    
    # Establish connection to SQLite database
    conn = sqlite3.connect(db_file)
    
    # Adjust the column name to match the hits CSV file and process dbcan_id values if needed
    hits_ids_column = f"{db_name}_id"
    if db_name == "dbcan":
        chunk[hits_ids_column] = chunk[hits_ids_column].str.replace(".hmm", "")
    ids = chunk[hits_ids_column].unique()
    
    # Construct the query based on the database name
    query = f"SELECT {ids_column}, {descriptions_column} FROM {table_name} WHERE {ids_column} IN ({','.join(['?'] * len(ids))})"
    if db_name == "dbcan":
        query = f"SELECT {ids_column}, {descriptions_column}, ec FROM {table_name} WHERE {ids_column} IN ({','.join(['?'] * len(ids))})"
    
    cursor = conn.cursor()
    cursor.execute(query, ids)
    results = cursor.fetchall()
    
    # Update chunk with descriptions and process according to the database type
    if db_name == "dbcan":
        # For dbcan, also fetch and format EC numbers
        descriptions_dict = {row[0]: (row[1], row[2]) for row in results}
        chunk[f"{db_name}_description"] = chunk[hits_ids_column].map(lambda x: descriptions_dict.get(x, ("", ""))[0])
        chunk["dbcan_EC"] = chunk[hits_ids_column].map(lambda x: format_dbcan_EC(descriptions_dict.get(x, ("", ""))[1]))
    elif db_name == "kegg":
        # For other databases, just update with descriptions
        descriptions_dict = {row[0]: row[1] for row in results}
        chunk[f"{db_name}_description"] = chunk[hits_ids_column].map(lambda x: descriptions_dict.get(x, ""))
        chunk["kegg_orthology"] = chunk[f"{db_name}_description"].apply(extract_kegg_orthology)
        chunk["kegg_EC"] = chunk[f"{db_name}_description"].apply(extract_kegg_EC)
        
        # Ensure the 'kegg_id' column exists if your logic does not create it earlier
        if 'kegg_id' not in chunk.columns:
            chunk['kegg_id'] = ''  # or any default logic you have for generating 'kegg_id'

        # Rename 'kegg_id' to 'kegg_gene_name'
        chunk.rename(columns={"kegg_id": "kegg_gene_name"}, inplace=True)
        
        # After extracting for kegg, rename "kegg_orthology" to "kegg_id"
        chunk.rename(columns={"kegg_orthology": "kegg_id"}, inplace=True)
    else:
        descriptions_dict = {row[0]: row[1] for row in results}
        chunk[f"{db_name}_description"] = chunk[hits_ids_column].map(lambda x: descriptions_dict.get(x, ""))

        
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

def format_dbcan_EC(ec_string):
    # Format EC numbers from the dbcan ec column with "EC:" prefix and semicolon separation
    if ec_string:
        ec_numbers = ec_string.split(',')  # Split the EC numbers string
        formatted_ec_numbers = '; '.join([f"EC:{ec}" for ec in ec_numbers if ec])  # Add "EC:" prefix and join with semicolon
        return formatted_ec_numbers
    else:
        return ""

def extract_kegg_EC(description):
    # Extract and format EC numbers from the description with "EC:" prefix and semicolon separation
    ec_start = description.find("[EC:")
    if ec_start != -1:
        ec_end = description.find("]", ec_start)
        ec_text = description[ec_start + 4:ec_end]  # Skip the "[EC:" part
        ec_numbers = re.findall(r'\b\d+\.\d+\.\d+\.\d+\b', ec_text)  # Extract EC numbers using regex
        formatted_ec_numbers = '; '.join([f"EC:{ec}" for ec in ec_numbers])  # Add "EC:" prefix and join with semicolon
        return formatted_ec_numbers
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