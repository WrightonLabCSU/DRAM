import json
import pandas as pd
import sys
from Bio import Phylo
import subprocess
import os

def run_guppy(jplace_path, output_dir):
    subprocess.run(['guppy', 'tog', jplace_path, '--output', f"{output_dir}/tree_with_placements.newick"], check=True)
    subprocess.run(['guppy', 'edpl', '--csv', jplace_path, '--output', f"{output_dir}/edpl.csv"], check=True)
    return f"{output_dir}/tree_with_placements.newick", f"{output_dir}/edpl.csv"

def load_phylogenetic_tree(tree_file):
    return Phylo.read(tree_file, 'newick')

def extract_closest_matches(jplace_path):
    with open(jplace_path, 'r') as file:
        data = json.load(file)
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = min(placement['p'], key=lambda x: x[3])  # Minimize by likelihood score
            edge_number, like_weight_ratio, likelihood = most_likely_placement[1], most_likely_placement[2], most_likely_placement[3]
            placements[gene_name] = (edge_number, like_weight_ratio, likelihood)
            print(f"Closest match for {gene_name}: Edge {edge_number}, LWR: {like_weight_ratio}, Likelihood: {likelihood}")
    return placements

def update_annotations(annotations_path, placements, tree, tree_mapping):
    annotations_df = pd.read_csv(annotations_path, sep='\t')
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

def main(jplace_file, mapping_file, annotations_file, output_file):
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    tree_path, edpl_path = run_guppy(jplace_path, output_dir)
    tree_mapping = load_tree_mapping(mapping_file)
    placements = extract_closest_matches(jplace_file)
    updated_annotations = update_annotations(annotations_file, placements, tree, tree_mapping)
    updated_annotations.to_csv(output_file, sep='\t', index=False)
    print("Updated annotations written to file.")

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python script.py <jplace_file> <mapping_file> <annotations_file> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
