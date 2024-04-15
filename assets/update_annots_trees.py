import json
import sys
import re

def load_jplace_file(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def parse_newick_tree(tree_string):
    node_stack = []
    parent_map = {}
    label_map = {}
    buf = ""
    node_id = 0  # Assign an ID to every node for tracking
    for char in tree_string:
        if char in "(,);":
            if buf:
                label_match = re.match(r'([^\[\]\:\{\}]+)\{(\d+)\}', buf)
                if label_match:
                    label, num = label_match.groups()
                    label_map[int(num)] = label
                node_stack.append((buf, node_id))  # Push node with its ID
                node_id += 1
                buf = ""
            if char == ")":  # Node is complete
                node, nid = node_stack.pop()
                if node_stack:
                    parent_map[nid] = node_stack[-1][1]  # Map to parent ID
            if char == "(":
                buf = char  # Include the character to handle nested trees
        else:
            buf += char
    print("Label Map:", label_map)  # Debugging output
    print("Parent Map:", parent_map)  # Debugging output
    return parent_map, label_map

def find_closest_labeled_ancestor(node, parent_map, label_map):
    while node in parent_map:
        node = parent_map[node]
        if node in label_map:
            return label_map[node]
    return "No labeled ancestor found"

def extract_placement_details(jplace_data, parent_map, label_map):
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
