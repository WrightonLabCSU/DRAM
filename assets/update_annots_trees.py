import json
import sys
from Bio import Phylo
from io import StringIO

def load_jplace_file(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def load_and_parse_tree(tree_data):
    handle = StringIO(tree_data)
    try:
        tree = Phylo.read(handle, "newick")
        print("Tree loaded successfully.")
    except Exception as e:
        print(f"Failed to load tree: {e}")
        tree = None
    return tree

def find_closest_labeled_ancestor(clade, tree):
    try:
        path = tree.get_path(clade)
        for ancestor in reversed(path):
            if ancestor.name:
                return ancestor.name
        return "No labeled ancestor found"
    except Exception as e:
        print(f"Error finding ancestor: {e}")
        return "Error in path finding"

def extract_placement_details(jplace_data, tree):
    placements = jplace_data['placements']
    results = []
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            try:
                clades = list(tree.find_clades({"comment": str(edge_num)}))
                if clades:
                    clade = clades[0]
                    closest_leaf = find_closest_labeled_ancestor(clade, tree)
                else:
                    closest_leaf = "Clade not found"
                    print(f"No clades found for edge number: {edge_num}")
            except Exception as e:
                closest_leaf = "Error in clade finding"
                print(f"Error finding clades for edge number {edge_num}: {e}")

            for name, _ in placement['nm']:
                results.append((name, edge_num, likelihood, closest_leaf))
    return results

def print_placements(results):
    if results:
        for name, edge_num, likelihood, leaf in results:
            print(f"Sequence: {name}, Edge: {edge_num}, Likelihood: {likelihood}, Closest Leaf: {leaf}")
    else:
        print("No results to display.")

def main():
    if len(sys.argv) != 2:
        print("Usage: python update_annots_trees.py <jplace_path>")
        return
    
    jplace_path = sys.argv[1]
    jplace_data = load_jplace_file(jplace_path)
    if 'tree' not in jplace_data:
        print("No tree data found in the jplace file.")
        return
    tree = load_and_parse_tree(jplace_data['tree'])
    if tree:
        results = extract_placement_details(jplace_data, tree)
        print_placements(results)
    else:
        print("Failed to process the tree data.")

if __name__ == "__main__":
    main()
