import json
import sys
import subprocess
import os
import re
from Bio import Phylo

def find_label_for_edge(tree, edge_number):
    for clade in tree.find_clades():
        if str(clade.comment) == str(edge_number):
            return clade.name, clade.branch_length
    return "No matching label found", None

def load_phylogenetic_tree(tree_file):
    try:
        tree = Phylo.read(tree_file, 'newick')
        return tree
    except FileNotFoundError as e:
        print(f"Error: Tree file '{tree_file}' not found.")
        raise
    except Phylo.NewickParserError as e:
        print(f"Error parsing tree file '{tree_file}': {e}")
        raise

def run_guppy(jplace_file, output_dir):
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Run guppy tog command to generate tree_with_placements.newick file
        tree_output_path = os.path.join(output_dir, "tree_with_placements.newick")
        subprocess.run(['guppy', 'tog', jplace_file, '-o', tree_output_path], check=True)
        
        # Run guppy edpl command to generate edpl.csv file
        edpl_output_path = os.path.join(output_dir, "edpl.csv")
        subprocess.run(['guppy', 'edpl', '--csv', jplace_file, '-o', edpl_output_path], check=True)
        
        return tree_output_path, edpl_output_path
    except subprocess.CalledProcessError as e:
        print(f"Error running guppy command: {e}")
        raise
    except OSError as e:
        print(f"Error creating output directory: {e}")
        raise

def extract_tree_and_placements(jplace_file):
    with open(jplace_file, 'r') as file:
        data = json.load(file)
    
    # Extract tree from .jplace file
    tree_path = data['tree']

    # Load the tree from the file
    tree = load_phylogenetic_tree(tree_path)

    # Extract placements
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            placements[gene_name] = placement['p'][0][1]  # Extract edge number
    
    return tree, placements


def find_closest_tip_labels(tree, placements):
    closest_tip_labels = {}
    for gene_id, edge_number in placements.items():
        closest_tip_label, edge_length = find_label_for_edge(tree, edge_number)
        closest_tip_labels[gene_id] = closest_tip_label
    return closest_tip_labels

def print_closest_tip_labels(closest_tip_labels):
    for gene_id, tip_label in closest_tip_labels.items():
        print(f"Sample: {gene_id}, Closest Tip Label: {tip_label}")

def main(jplace_file):
    output_dir = './output'
    tree_path, _ = run_guppy(jplace_file, output_dir)
    tree, placements = extract_tree_and_placements(jplace_file)
    closest_tip_labels = find_closest_tip_labels(tree, placements)
    print_closest_tip_labels(closest_tip_labels)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python script.py <jplace_file>")
        sys.exit(1)
    main(sys.argv[1])
