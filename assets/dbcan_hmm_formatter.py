import pandas as pd
import argparse

def remove_extension(target_id):
    # Remove the ".hmm" extension from target_id
    return target_id.replace(".hmm", "")

def read_ch_dbcan_fam(file_path):
    # Read ch_dbcan_fam file with flexible delimiter handling
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line based on tabs
            fields = line.strip().split('\t', 1)
            lines.append(fields)
    return pd.DataFrame(lines, columns=['AA', 'Description'])

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="DBCAN HMM Formatter")
    parser.add_argument("--hits_csv", required=True, help="Path to hits CSV file")
    parser.add_argument("--fam", required=True, help="Path to ch_dbcan_fam file")
    parser.add_argument("--subfam", required=True, help="Path to ch_dbcan_subfam file")
    parser.add_argument("--output", required=True, help="Output file name")

    args = parser.parse_args()

    # Read hits file
    hits_df = pd.read_csv(args.hits_csv, sep=",")  # Assuming comma as the delimiter

    # Debugging: Print column names and contents of hits_df
    print("Column names of hits_df:", hits_df.columns)
    print("Contents of hits_df:")
    print(hits_df.head())


    # Read ch_dbcan_fam file
    ch_dbcan_fam = read_ch_dbcan_fam(args.fam)

    # Read ch_dbcan_subfam file
    ch_dbcan_subfam = pd.read_csv(args.subfam, sep="\t", header=None, names=['AA', 'GenBank', 'EC'])

    # Remove ".hmm" extension from target_id in hits_df
    hits_df['target_id'] = hits_df['target_id'].apply(remove_extension)

    # Merge hits_df with ch_dbcan_fam to get subfamily
    hits_df = pd.merge(hits_df, ch_dbcan_fam, left_on='target_id', right_on='AA', how='left')

    # Merge hits_df with ch_dbcan_subfam to get GenBank and EC
    merged_df = pd.merge(hits_df, ch_dbcan_subfam, on='AA', how='left')

    # Group by query_id and aggregate multiple GenBank and EC values
    grouped_df = merged_df.groupby('query_id').agg({
        'target_id': 'first',
        'bitScore': 'first',
        'subfamily': 'first',
        'GenBank': lambda x: "; ".join(x.dropna()),  # Concatenate GenBank values
        'EC': lambda x: "; ".join(x.dropna())  # Concatenate EC values
    }).reset_index()

    # Save the formatted hits to the output file
    grouped_df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
