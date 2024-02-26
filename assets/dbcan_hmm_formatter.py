import pandas as pd
import argparse
import sqlite3

def calculate_rank(row):
    return row['score_rank'] if 'score_rank' in row and row['full_score'] > row['score_rank'] else row['full_score']

def calculate_bit_score(row):
    return row['full_score'] / row['domain_number']

def calculate_coverage(row):
    return (row['target_end'] - row['target_start']) / row['target_length']

def fetch_descriptions_from_db(target_ids, db_file):
    conn = sqlite3.connect(db_file)
    descriptions = {}
    for target_id in target_ids:
        cursor = conn.execute("SELECT description, ec FROM dbcan_description WHERE id=?", (target_id,))
        row = cursor.fetchone()
        if row:
            description, ec = row
            descriptions[target_id] = {'description': description, 'ec': ec}
        else:
            descriptions[target_id] = {'description': "", 'ec': ""}  # Handle case where description is not found
    conn.close()
    return descriptions

def main():
    parser = argparse.ArgumentParser(description="Format HMM search results.")
    parser.add_argument("--hits_csv", type=str, help="Path to the HMM search results CSV file.")
    parser.add_argument("--output", type=str, help="Path to the formatted output file.")
    parser.add_argument("--db_file", type=str, help="Path to the SQLite database file.")

    args = parser.parse_args()

    print("Loading HMM search results CSV file...")
    hits_df = pd.read_csv(args.hits_csv)
    print(f"Loaded HMM search results from: {args.hits_csv}")

    print("First few lines of hits_df:")
    print(hits_df.head())

    required_columns = ['query_id', 'query_start', 'query_end', 'strandedness', 'target_id', 'score_rank', 'full_score', 'domain_number', 'target_length', 'target_start', 'target_end']
    missing_columns = [col for col in required_columns if col not in hits_df.columns]
    
    if missing_columns:
        print(f"Error: Missing columns in hits_df: {missing_columns}")
        return

    print("Processing HMM search results...")
    hits_df['target_id'] = hits_df['target_id'].str.replace(r'.hmm', '', regex=True)

    hits_df['bitScore'] = hits_df.apply(calculate_bit_score, axis=1)
    print("Bit scores calculated.")

    hits_df['score_rank'] = hits_df.apply(calculate_rank, axis=1)
    print("Ranks calculated.")

    # Calculate coverage
    hits_df['perc_cov'] = hits_df.apply(calculate_coverage, axis=1)
    print("Coverage calculated.")

    hits_df.dropna(subset=['score_rank'], inplace=True)
    print("NaN values dropped.")

    # Save the intermediate DataFrame with required columns for step 1
    intermediate_columns = ['query_id', 'query_start', 'query_end', 'strandedness', 'target_id', 'score_rank', 'bitScore']
    intermediate_df = hits_df[intermediate_columns].copy()
    print("Intermediate DataFrame created with required columns.")

    # Save the intermediate DataFrame to a CSV file for step 1
    intermediate_output = args.output.split('.')[0] + "_intermediate.csv"
    try:
        intermediate_df.to_csv(intermediate_output, index=False)
        print(f"Intermediate output saved to: {intermediate_output}")
    except Exception as e:
        print(f"Error occurred while saving the intermediate output: {e}")
        return

    print("Intermediate process completed successfully.")

    # Step 2: Add descriptions and ECs to the intermediate DataFrame
    print("Fetching descriptions and ECs from the database...")
    target_ids = intermediate_df['target_id'].unique()
    descriptions = fetch_descriptions_from_db(target_ids, args.db_file)

    # Assign descriptions and ECs to hits
    intermediate_df['dbcan_description'] = intermediate_df['target_id'].map(lambda x: descriptions[x]['description'])
    intermediate_df['dbcan_ec'] = intermediate_df['target_id'].map(lambda x: descriptions[x]['ec'])
    print("Descriptions and ECs assigned to hits.")

    # Save the final formatted output to CSV
    selected_columns = ['query_id', 'query_start', 'query_end', 'strandedness', 'dbcan_id', 'dbcan_score_rank', 'dbcan_bitScore', 'dbcan_description', 'dbcan_ec']
    final_output_df = intermediate_df[selected_columns].copy()

    try:
        final_output_df.to_csv(args.output, index=False)
        print(f"Formatted output saved to: {args.output}")
    except Exception as e:
        print(f"Error occurred while saving the formatted output: {e}")

    print("Process completed successfully!")

if __name__ == "__main__":
    main()
