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
    return closest_leaf if closest_leaf else "No labeled ancestor found"

def extract_placement_details(jplace_data, tree):
    placements = jplace_data['placements']
    results = {}
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            clades = list(tree.find_clades({"comment": f"&edge={edge_num}"}))
            closest_leaf = "Clade not found" if not clades else find_closest_labeled_ancestor(clades[0], tree)
            for name, _ in placement['nm']:
                results[name] = closest_leaf
    return results

def update_tsv_file(tsv_path, placements, output_path):
    df = pd.read_csv(tsv_path, sep='\t')
    df['tree_verified'] = df['query_id'].map(placements).fillna('No labeled ancestor found')
    gene_number_idx = df.columns.get_loc('gene_number')
    columns = list(df.columns)
    columns.insert(gene_number_idx + 1, 'tree_verified')
    df = df[columns]
    df.to_csv(output_path, sep='\t', index=False)

def main():
    if len(sys.argv) < 4:
        print("Usage: python update_annots_trees.py <jplace_path> <tsv_path> <output_path>")
        return
    jplace_path, tsv_path, output_path = sys.argv[1:4]
    jplace_data = load_jplace_file(jplace_path)
    tree = load_and_parse_tree(jplace_data['tree'])
    placements = extract_placement_details(jplace_data, tree)
    update_tsv_file(tsv_path, placements, output_path)
    print("TSV file updated and saved to:", output_path)

if __name__ == "__main__":
    main()
