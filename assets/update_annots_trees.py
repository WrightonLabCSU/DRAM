import json
import pandas as pd
import sys
from Bio import Phylo
import subprocess
import os

def run_guppy(jplace_path, output_dir):
    # Generate tree with placements
    tree_path = f"{output_dir}/tree_with_placements.newick"
    subprocess.run(['guppy', 'tog', '--xml', jplace_path, '-o', tree_path], check=True)
    print(f"Tree with placements written to {tree_path}")

    # Calculate EDPL
    edpl_path = f"{output_dir}/edpl.csv"
    subprocess.run(['guppy', 'edpl', '--point-mass', '--csv', jplace_path, '-o', edpl_path], check=True)
    print(f"EDPL values written to {edpl_path}")

    return tree_path, edpl_path

def load_phylogenetic_tree(tree_file):
    """Load the phylogenetic tree from a Newick file."""
    return Phylo.read(tree_file, 'newick')

def extract_placements(jplace_path):
    """Extract placement edge numbers from the jplace file."""
    with open(jplace_path, 'r') as file:
        data = json.load(file)
    
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = placement['p'][0]
            edge_number = most_likely_placement[1]
            placements[gene_name] = edge_number
            print(f"Extracted placement for {gene_name}: Edge {edge_number}")
    return placements

def load_tree_mapping(mapping_path):
    mapping_df = pd.read_csv(mapping_path, sep='\t')
    tree_map = mapping_df.set_index('gene')['call'].to_dict()
    print(f"Loaded tree mapping for {len(tree_map)} genes.")
    return tree_map

def find_named_ancestor(tree, edge_number, tree_mapping):
    """Find the nearest named ancestor of a given edge in the tree."""
    try:
        node = next(tree.find_clades({"name": str(edge_number)}))
    except StopIteration:
        print(f"No node found for edge number {edge_number}")
        return None

    # Trace back to find a named ancestor
    while node:
        if node.name and node.name in tree_mapping:
            print(f"Found matching ancestor: {node.name} for edge {edge_number}")
            return tree_mapping[node.name]
        node = node.parent
        if node:
            print(f"Visiting node: {node.name}")

    print(f"No matching ancestor found in mapping for edge number {edge_number}")
    return None

def update_annotations(annotations_path, placements, tree, tree_mapping):
    annotations_df = pd.read_csv(annotations_path, sep='\t')
    annotations_df['tree-verified'] = None  # Initialize column
    print(f"Attempting to update {len(annotations_df)} annotations from placements.")

    for index, row in annotations_df.iterrows():
        gene_id = row['query_id']
        if gene_id in placements:
            edge = placements[gene_id]
            ancestor = find_named_ancestor(tree, edge, tree_mapping)
            if ancestor:
                annotations_df.at[index, 'tree-verified'] = ancestor
                print(f"Gene {gene_id} placed on edge {edge} which corresponds to {ancestor}")
            else:
                print(f"Gene {gene_id} placed on edge {edge} but no matching tree mapping found.")

    return annotations_df

def main(jplace_file, mapping_file, annotations_file, output_file, tree_file):
    # Directory to store output files
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    
    # Run guppy to get the tree and EDPL
    tree_path, edpl_path = run_guppy(jplace_file, output_dir)
    
    # Load the phylogenetic tree with placements
    tree = Phylo.read(tree_path, 'newick')
    
    # Load tree mappings and placements
    tree_mapping = load_tree_mapping(mapping_file)
    placements = extract_placements(jplace_file)
    
    # Update annotations
    updated_annotations = update_annotations(annotations_file, placements, tree, tree_mapping)
    updated_annotations.to_csv(output_file, sep='\t', index=False)
    print("Updated annotations written to file.")

if __name__ == '__main__':
    jplace_file = sys.argv[1]
    mapping_file = sys.argv[2]
    annotations_file = sys.argv[3]
    output_file = sys.argv[4]
    tree_file = sys.argv[5]
    main(jplace_file, mapping_file, annotations_file, output_file, tree_file)
