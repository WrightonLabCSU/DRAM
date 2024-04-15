import json
import sys
import re
import sqlite3
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

def extract_placement_details(jplace_data, tree):
    placements = jplace_data['placements']
    placement_map = {}
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            clades = list(tree.find_clades({"comment": f"&edge={edge_num}"}))
            if clades:
                clade = clades[0]
                closest_leaf = find_closest_labeled_ancestor(clade, tree)
            else:
                closest_leaf = "No labeled ancestor found"
            for name, _ in placement['nm']:
                placement_map[name] = closest_leaf if closest_leaf != "No labeled ancestor found" else ""
    return placement_map

def update_database_and_tsv(tsv_path, db_path, output_db_path, output_tsv_path, placement_map):
    # Read TSV into DataFrame
    df = pd.read_csv(tsv_path, sep='\t')
    # Map placements to DataFrame
    df['tree_verified'] = df['query_id'].map(placement_map)
    
    # Write to new TSV
    df.to_csv(output_tsv_path, sep='\t', index=False)
    
    # Update SQLite Database
    conn = sqlite3.connect(db_path)
    df.to_sql('annotations', conn, if_exists='replace', index=False)
    conn.close()

    # Export updated DataFrame to new SQLite DB
    conn_out = sqlite3.connect(output_db_path)
    df.to_sql('annotations', conn_out, if_exists='replace', index=False)
    conn_out.close()

def main():
    if len(sys.argv) != 6:
        print("Usage: python update_annots_trees.py <jplace_path> <tsv_path> <db_path> <output_db_path> <output_tsv_path>")
        return
    
    jplace_path, tsv_path, db_path, output_db_path, output_tsv_path = sys.argv[1:]
    jplace_data = load_jplace_file(jplace_path)
    tree = load_and_parse_tree(jplace_data['tree'])
    placement_map = extract_placement_details(jplace_data, tree)
    update_database_and_tsv(tsv_path, db_path, output_db_path, output_tsv_path, placement_map)
    print("Updates complete. Data saved to", output_db_path, "and", output_tsv_path)

if __name__ == "__main__":
    main()
