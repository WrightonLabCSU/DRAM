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
                    closest_leaf = f"{tree_mapping[closest_leaf]};{closest_leaf}"
                else:
                    closest_leaf = f"No mapping found;{closest_leaf}"
            else:
                closest_leaf = f"Clade not found;{edge_num}"
            for name, _ in placement['nm']:
                placement_map[name] = closest_leaf
    return placement_map

def update_tsv(tsv_path, output_tsv_path, placement_map):
    df = pd.read_csv(tsv_path, sep='\t')
    df['tree_verified'] = df['query_id'].map(placement_map).fillna('')

    # Reorder columns to place 'tree_verified' after 'gene_number'
    col_order_start = df.columns.tolist()[:df.columns.get_loc('gene_number')+1] + ['tree_verified']
    col_order_end = [col for col in df.columns if col not in col_order_start]
    df = df[col_order_start + col_order_end]

    df.to_csv(output_tsv_path, sep='\t', index=False)

def main():
    if len(sys.argv) != 5:
        print("Usage: python update_annots_trees.py <jplace_path> <tsv_path> <mapping_tsv> <output_tsv_path>")
        sys.exit(1)

    jplace_path, tsv_path, mapping_tsv, output_tsv_path = sys.argv[1:]
    jplace_data = load_jplace_file(jplace_path)
    tree = load_and_parse_tree(jplace_data['tree'])
    tree_mapping = load_tree_mapping(mapping_tsv)
    placement_map = extract_placement_details(jplace_data, tree, tree_mapping)
    update_tsv(tsv_path, output_tsv_path, placement_map)

if __name__ == "__main__":
    main()
