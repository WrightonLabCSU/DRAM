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

def extract_placements(jplace_file):
    with open(jplace_file, 'r') as file:
        data = json.load(file)
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = min(placement['p'], key=lambda x: x[3])  # Minimize by likelihood score
            edge_number, _, _ = most_likely_placement[1], most_likely_placement[2], most_likely_placement[3]
            placements[gene_name] = edge_number
    return placements

def load_tree_mapping(mapping_file):
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    return mapping_df.set_index('gene')['call'].to_dict()

def find_closest_tip_label(tree, placements):
    closest_tip_labels = {}
    for gene_id, edge_number in placements.items():
        closest_tip_labels[gene_id] = find_label_for_edge(tree, edge_number)
    return closest_tip_labels

def find_label_for_edge(tree, edge_number):
    for clade in tree.find_clades():
        if str(clade.comment) == str(edge_number):
            return clade.name
    return "No matching label found"

def print_closest_tip_labels(closest_tip_labels):
    for gene_id, tip_label in closest_tip_labels.items():
        print(f"Sample: {gene_id}, Closest Tip Label: {tip_label}")

def main(jplace_file, mapping_file):
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    tree_path, _ = run_guppy(jplace_file, output_dir)
    tree = load_phylogenetic_tree(tree_path)
    placements = extract_placements(jplace_file)
    tree_mapping = load_tree_mapping(mapping_file)
    closest_tip_labels = find_closest_tip_label(tree, placements)
    print_closest_tip_labels(closest_tip_labels)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python script.py <jplace_file> <mapping_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
