import json
import sys
from Bio import Phylo

def load_jplace_file(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def load_and_parse_tree(tree_data):
    from io import StringIO
    handle = StringIO(tree_data)
    tree = Phylo.read(handle, "newick")
    return tree

def find_closest_labeled_ancestor(clade, tree):
    # Find the path from this clade to the root
    path = tree.get_path(clade)
    # Traverse from clade to root, finding the first labeled ancestor
    for ancestor in reversed(path):
        if ancestor.name:
            return ancestor.name
    return "No labeled ancestor found"

def extract_placement_details(jplace_data, tree):
    placements = jplace_data['placements']
    results = []
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            clade = tree.find_clades({"comment": str(edge_num)}).next()  # Assumes edge_num is stored in comment
            closest_leaf = find_closest_labeled_ancestor(clade, tree)
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
