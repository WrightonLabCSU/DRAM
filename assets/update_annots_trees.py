import json
import pandas as pd
import sys
from Bio import Phylo
import subprocess
import os

def run_guppy(jplace_path, output_dir):
    tree_path = f"{output_dir}/tree_with_placements.newick"
    subprocess.run(['guppy', 'tog', jplace_path, '-o', tree_path], check=True)
    print(f"Tree with placements written to {tree_path}")

    edpl_path = f"{output_dir}/edpl.csv"
    subprocess.run(['guppy', 'edpl', '--csv', jplace_path, '-o', edpl_path], check=True)
    print(f"EDPL values written to {edpl_path}")

    return tree_path, edpl_path

def load_phylogenetic_tree(tree_file):
    return Phylo.read(tree_file, 'newick')

def extract_placements(jplace_path):
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
    target_node = next((node for node in tree.find_clades() if str(node.name) == str(edge_number)), None)
    if target_node is None:
        print(f"No node found for edge number {edge_number}")
        return None

    # Attempt to find a named ancestor directly in the path to the root
    ancestor = trace_to_root_for_named_ancestor(target_node, tree_mapping)
    if ancestor:
        return ancestor

    # If no named ancestor found in the path to the root, check nearby nodes
    nearest_node, distance = find_nearest_node(tree, target_node)
    if nearest_node:
        print(f"Nearest node to {target_node.name} is {nearest_node.name} at distance {distance}")
        return trace_to_root_for_named_ancestor(nearest_node, tree_mapping)

    return None

def trace_to_root_for_named_ancestor(node, tree_mapping):
    while node:
        if node.name and node.name in tree_mapping:
            print(f"Found matching ancestor: {node.name} for edge {node.name}")
            return tree_mapping[node.name]
        node = node.parent
    return None


def find_nearest_node(tree, target_node):
    min_distance = float('inf')
    nearest_node = None
    for node in tree.get_terminals():
        if node != target_node:
            distance = tree.distance(target_node, node)
            if distance < min_distance:
                min_distance = distance
                nearest_node = node
    return nearest_node, min_distance

def update_annotations(annotations_path, placements, tree, tree_mapping):
    annotations_df = pd.read_csv(annotations_path, sep='\t')
    annotations_df['tree-verified'] = None
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
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    tree_path, edpl_path = run_guppy(jplace_file, output_dir)
    tree = load_phylogenetic_tree(tree_path)
    tree_mapping = load_tree_mapping(mapping_file)
    placements = extract_placements(jplace_file)
    updated_annotations = update_annotations(annotations_file, placements, tree, tree_mapping)
    updated_annotations.to_csv(output_file, sep='\t', index=False)
    print("Updated annotations written to file.")

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("Usage: python script.py <jplace_file> <mapping_file> <annotations_file> <output_file> <tree_file>")
        sys.exit(1)
    jplace_file = sys.argv[1]
    mapping_file = sys.argv[2]
    annotations_file = sys.argv[3]
    output_file = sys.argv[4]
    tree_file = sys.argv[5]
    main(jplace_file, mapping_file, annotations_file, output_file, tree_file)
