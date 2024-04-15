import json
import pandas as pd
import sys
from Bio import Phylo
import subprocess
import os

def run_guppy(jplace_file, output_dir):
    subprocess.run(['guppy', 'tog', jplace_file, '-o', f"{output_dir}/tree_with_placements.newick"], check=True)
    subprocess.run(['guppy', 'edpl', '--csv', jplace_file, '-o', f"{output_dir}/edpl.csv"], check=True)
    return f"{output_dir}/tree_with_placements.newick", f"{output_dir}/edpl.csv"

def load_phylogenetic_tree(tree_file):
    return Phylo.read(tree_file, 'newick')

def extract_closest_matches(jplace_file, tree):
    with open(jplace_file, 'r') as file:
        data = json.load(file)
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = min(placement['p'], key=lambda x: x[3])  # Minimize by likelihood score
            edge_number, like_weight_ratio, likelihood = most_likely_placement[1], most_likely_placement[2], most_likely_placement[3]
            closest_existing_label = find_closest_existing_label(tree, edge_number)
            placements[gene_name] = (edge_number, like_weight_ratio, likelihood, closest_existing_label)
            print(f"Closest match for {gene_name}: Edge {edge_number}, LWR: {like_weight_ratio}, Likelihood: {likelihood}, Closest label: {closest_existing_label}")
    return placements

def find_closest_existing_label(tree, edge_number):
    # This function needs to be refined to effectively find the closest existing tip or node.
    for clade in tree.find_clades():
        if clade.name and str(clade.name) == str(edge_number):
            return clade.name
    return "No matching label found"

def main(jplace_file, output_file):
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    tree_path, edpl_path = run_guppy(jplace_file, output_dir)
    tree = load_phylogenetic_tree(tree_path)
    placements = extract_closest_matches(jplace_file, tree)
    # This will print detailed information for debugging

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python script.py <jplace_file> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

    
def update_annotations(annotations_file, placements, tree_mapping):
    try:
        annotations_df = pd.read_csv(annotations_file, sep='\t')
    except Exception as e:
        print(f"Error reading annotations file: {e}")
        sys.exit(1)

    annotations_df['tree-verified'] = None
    for index, row in annotations_df.iterrows():
        gene_id = row['query_id']
        if gene_id in placements:
            edge, lwr, likelihood = placements[gene_id]
            if edge in tree_mapping:
                annotations_df.at[index, 'tree-verified'] = tree_mapping[edge]
                print(f"Gene {gene_id} closest match on edge {edge} mapped to {tree_mapping[edge]}")
            else:
                print(f"Gene {gene_id} placed on edge {edge} but no matching tree mapping found.")
    return annotations_df
