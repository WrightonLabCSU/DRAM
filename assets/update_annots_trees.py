import json
import pandas as pd
import sys
import subprocess
import os
from Bio import Phylo

def run_guppy(jplace_file, output_dir):
    try:
        subprocess.run(['guppy', 'tog', jplace_file, '-o', f"{output_dir}/tree_with_placements.newick"], check=True)
        subprocess.run(['guppy', 'edpl', '--csv', jplace_file, '-o', f"{output_dir}/edpl.csv"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running guppy: {e}")
        sys.exit(1)
    return f"{output_dir}/tree_with_placements.newick", f"{output_dir}/edpl.csv"

def load_phylogenetic_tree(tree_file):
    try:
        tree = Phylo.read(tree_file, 'newick')
    except Exception as e:
        print(f"Error reading tree file: {e}")
        sys.exit(1)
    return tree

def extract_closest_matches(jplace_file):
    try:
        with open(jplace_file, 'r') as file:
            data = json.load(file)
    except Exception as e:
        print(f"Error loading jplace file: {e}")
        sys.exit(1)

    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = min(placement['p'], key=lambda x: x[3])
            edge_number, like_weight_ratio, likelihood = most_likely_placement[1], most_likely_placement[2], most_likely_placement[3]
            placements[gene_name] = (edge_number, like_weight_ratio, likelihood)
            print(f"Closest match for {gene_name}: Edge {edge_number}, LWR: {like_weight_ratio}, Likelihood: {likelihood}")
    return placements

def load_tree_mapping(mapping_file):
    try:
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        tree_map = mapping_df.set_index('gene')['call'].to_dict()
        print(f"Loaded tree mapping for {len(tree_map)} genes.")
    except Exception as e:
        print(f"Error loading mapping file: {e}")
        sys.exit(1)
    return tree_map

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

def main(jplace_file, mapping_file, annotations_file, output_file):
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    tree_path, edpl_path = run_guppy(jplace_file, output_dir)
    tree_mapping = load_tree_mapping(mapping_file)
    placements = extract_closest_matches(jplace_file)
    updated_annotations = update_annotations(annotations_file, placements, tree_mapping)
    updated_annotations.to_csv(output_file, sep='\t', index=False)
    print("Updated annotations written to file.")

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python script.py <jplace_file> <mapping_file> <annotations_file> <output_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
