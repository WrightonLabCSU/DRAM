import json
import sys
import re
import pandas as pd
from Bio import Phylo
from io import StringIO

def load_jplace_file(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def load_and_parse_tree(tree_data):
    tree_data = re.sub(r'\{(\d+)\}', r'[&edge=\g<1>]', tree_data)
    handle = StringIO(tree_data)
    tree = Phylo.read(handle, "newick")
    return tree

def find_closest_labeled_ancestor(clade, tree):
    if clade.is_terminal() and clade.name:
        return clade.name
    path = tree.get_path(clade)
    for ancestor in reversed(path):
        if ancestor.is_terminal() and ancestor.name:
            return ancestor.name
    min_distance = float('inf')
    closest_leaf = None
    for leaf in tree.get_terminals():
        distance = tree.distance(clade, leaf)
        if distance < min_distance and leaf.name:
            min_distance = distance
            closest_leaf = leaf.name
    return closest_leaf if closest_leaf else ""

def load_tree_mapping(mapping_tsv):
    # Load mapping TSV into a dictionary
    df = pd.read_csv(mapping_tsv, sep='\t')
    return dict(zip(df['gene'], df['call']))

def extract_placement_details(jplace_data, tree, tree_mapping):
    placements = jplace_data['placements']
    placement_map = {}
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            clades = list(tree.find_clades({"comment": f"&edge={edge_num}"}))
            if clades:
                clade = clades[0]
                closest_leaf = find_closest_labeled_ancestor(clade, tree)
                if closest_leaf and closest_leaf in tree_mapping:
                    # Combine tree verified label with additional metadata
                    closest_leaf = f"{tree_mapping[closest_leaf]};{closest_leaf}"
            else:
                closest_leaf = ""
                print(f"No clades found for edge number: {edge_num}")
            for name, _ in placement['nm']:
                placement_map[name] = closest_leaf
    return placement_map

def update_tsv(tsv_path, output_tsv_path, placement_map):
    # Read TSV into DataFrame
    df = pd.read_csv(tsv_path, sep='\t')

    # Map placements to DataFrame
    df['tree_verified'] = df['query_id'].map(placement_map).fillna('')

    # Define the new order of columns, ensuring 'tree_verified' is after 'gene_number'
    new_column_order = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank', 'gene_number', 'tree_verified']
    
    # If there are additional columns, add them to the list in their original order
    additional_columns = [col for col in df.columns if col not in new_column_order]
    new_column_order.extend(additional_columns)
    
    # Reorder DataFrame according to new column order
    df = df[new_column_order]

    # Write to new TSV
    df.to_csv(output_tsv_path, sep='\t', index=False)

def main():
    if len(sys.argv) != 6:
        print("Usage: python update_annots_trees.py <jplace_path> <tsv_path> <mapping_tsv> <output_tsv_path>")
        return
    
    jplace_path, tsv_path, mapping_tsv, output_tsv_path = sys.argv[1:]
    jplace_data = load_jplace_file(jplace_path)
    tree = load_and_parse_tree(jplace_data['tree'])
    tree_mapping = load_tree_mapping(mapping_tsv)
    placement_map = extract_placement_details(jplace_data, tree, tree_mapping)
    update_tsv(tsv_path, output_tsv_path, placement_map)
    print("Updates complete. Data saved to", output_tsv_path)

if __name__ == "__main__":
    main()
