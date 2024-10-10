import pandas as pd
import argparse

def find_rbh(forward_file, reverse_file, output_file):
    # Load the forward and reverse search results
    forward_df = pd.read_csv(forward_file, sep='\t', header=None)
    reverse_df = pd.read_csv(reverse_file, sep='\t', header=None)
    
    # Assuming the columns are [query, target, bitscore,...], adjust if your format is different
    # Map reverse hits for easy lookup
    reverse_hits = {row[1]: row[0] for index, row in reverse_df.iterrows()}

    # Filter for reciprocal best hits
    rbh = []
    for index, row in forward_df.iterrows():
        query, target = row[0], row[1]
        # Check if the target of the forward hit is the query of a reverse hit, and it maps back to the original query
        if reverse_hits.get(query) == target:
            rbh.append(row)

    # Convert RBH list to DataFrame
    rbh_df = pd.DataFrame(rbh)
    if not rbh_df.empty:
        # Save the RBH to an output file
        rbh_df.to_csv(output_file, sep='\t', header=False, index=False)
    else:
        print("No reciprocal best hits found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter Reciprocal Best Hits (RBH) from forward and reverse MMseqs2 search results.")
    parser.add_argument('--forward', type=str, required=True, help='File path to forward search results.')
    parser.add_argument('--reverse', type=str, required=True, help='File path to reverse search results.')
    parser.add_argument('--output', type=str, required=True, help='Output file path for combined RBH results.')

    args = parser.parse_args()
    find_rbh(args.forward, args.reverse, args.output)
