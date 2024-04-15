import json
import sys
import re

def load_jplace_file(jplace_path):
    """ Load and parse the .jplace JSON file. """
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def parse_newick_tree(tree_string):
    """ Parse the Newick tree string to create a mapping of edge numbers to labels. """
    pattern = r'([^(,;]+)\{(\d+)\}'
    matches = re.findall(pattern, tree_string)
    return {int(num): label for label, num in matches}

def extract_placement_details(jplace_data, edge_to_label):
    """ Extract and print details from placements, including closest leaves. """
    placements = jplace_data['placements']
    results = []
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            closest_leaf = edge_to_label.get(edge_num, "Unknown leaf")
            for name, _ in placement['nm']:
                results.append((name, edge_num, likelihood, closest_leaf))
    return results

def print_placements(results):
    """ Print the placement details for each sequence. """
    for name, edge_num, likelihood, leaf in results:
        print(f"Sequence: {name}, Edge: {edge_num}, Likelihood: {likelihood}, Closest Leaf: {leaf}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python update_annots_trees.py <jplace_path>")
        return
    
    jplace_path = sys.argv[1]
    jplace_data = load_jplace_file(jplace_path)
    edge_to_label = parse_newick_tree(jplace_data['tree'])
    results = extract_placement_details(jplace_data, edge_to_label)
    print_placements(results)

if __name__ == "__main__":
    main()
