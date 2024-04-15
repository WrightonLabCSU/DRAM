import json
import sys
import re
from Bio import Phylo
from io import StringIO

def load_jplace_file(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def load_and_parse_tree(tree_data):
    # Convert edge numbers to comments within square brackets
    tree_data = re.sub(r'\{(\d+)\}', r'[&edge=\g<1>]', tree_data)
    handle = StringIO(tree_data)
    tree = Phylo.read(handle, "newick")
    return tree

def find_closest_labeled_ancestor(clade, tree):
    # Start with the given clade and check if it's a leaf and labeled
    if clade.is_terminal() and clade.name:
        return clade.name

    # If not a leaf, get the path to the root and look for the first labeled leaf in the path
    path = tree.get_path(clade)
    for ancestor in reversed(path):  # Start checking from the clade upwards
        if ancestor.is_terminal() and ancestor.name:
            return ancestor.name

    # As a fallback, search all leaves and find the closest one if no labeled ancestors found in the path
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
    results = []
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            clades = list(tree.find_clades({"comment": f"&edge={edge_num}"}))
            if clades:
                clade = clades[0]
                closest_leaf = find_closest_labeled_ancestor(clade, tree)
            else:
                closest_leaf = "Clade not found"
                print(f"No clades found for edge number: {edge_num}")
            for name, _ in placement['nm']:
                results.append((name, edge_num, likelihood, closest_leaf))
    return results

def print_placements(results):
    for name, edge_num, likelihood, leaf in results:
        print(f"Sequence: {name}, Edge: {edge_num}, Likelihood: {likelihood}, Closest Leaf: {leaf}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python update_annots_trees.py <jplace_path>")
        return

    jplace_path = sys.argv[1]
    jplace_data = load_jplace_file(jplace_path)
    tree = load_and_parse_tree(jplace_data['tree'])
    results = extract_placement_details(jplace_data, tree)
    print_placements(results)

if __name__ == "__main__":
    main()
