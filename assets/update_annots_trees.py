import json
import sys
import re

def load_jplace_file(jplace_path):
    """ Load and parse the .jplace JSON file. """
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def parse_newick_tree(tree_string):
    """ Parse the Newick tree to create mappings of nodes to their parents and labels. """
    node_stack = []
    parent_map = {}
    label_map = {}
    buf = ""
    for char in tree_string:
        if char in "(,);":
            if buf:
                label_match = re.match(r'([^\[\]\:\{\}]+)\{(\d+)\}', buf)
                if label_match:
                    label, node_id = label_match.groups()
                    label_map[int(node_id)] = label
                node_stack.append(buf)
                buf = ""
            if char == ")":  # End of a subtree
                node = node_stack.pop()
                if node_stack:
                    parent_map[node] = node_stack[-1]
            if char in "(,;":  # Start of a new node or subtree
                if node_stack:
                    parent_map[buf] = node_stack[-1]
        else:
            buf += char
    return parent_map, label_map

def find_closest_labeled_ancestor(node, parent_map, label_map):
    """ Traverse up from a node to find the closest ancestor with a label. """
    while node in parent_map:
        node = parent_map[node]
        if node in label_map:
            return label_map[node]
    return "No labeled ancestor found"

def extract_placement_details(jplace_data, parent_map, label_map):
    """ Extract and print details from placements, including closest leaves. """
    placements = jplace_data['placements']
    results = []
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            closest_leaf = label_map.get(edge_num, find_closest_labeled_ancestor(edge_num, parent_map, label_map))
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
    parent_map, label_map = parse_newick_tree(jplace_data['tree'])
    results = extract_placement_details(jplace_data, parent_map, label_map)
    print_placements(results)

if __name__ == "__main__":
    main()
