import json
import sys
import subprocess
import re

def find_label_for_edge(tree, edge_number):
    # Parse the Newick format tree string
    clade_regex = re.compile(r'(\([^()]+)\)?(:[-+]?[0-9]*\.?[0-9]+)?([,)]{1})')
    node_regex = re.compile(r'([^():,]+)(:.+)?')
    nodes = clade_regex.findall(tree)
    
    # Traverse the tree to find the label associated with the given edge number
    for node in nodes:
        matches = node_regex.findall(node[0])
        label = matches[0][0]
        if len(matches[0]) > 1:  # Check if node has edge length
            edge_length = matches[0][1][1:]  # Remove ':' from edge length
            if edge_length:
                edge_length = float(edge_length)
        else:
            edge_length = None
        
        if edge_length == edge_number:
            return label
    
    return "No matching label found"
def run_guppy(jplace_file, output_dir):
    subprocess.run(['guppy', 'tog', jplace_file, '-o', f"{output_dir}/tree_with_placements.newick"], check=True)
    subprocess.run(['guppy', 'edpl', '--csv', jplace_file, '-o', f"{output_dir}/edpl.csv"], check=True)
    return f"{output_dir}/tree_with_placements.newick", f"{output_dir}/edpl.csv"

def extract_tree_and_placements(jplace_file):
    with open(jplace_file, 'r') as file:
        data = json.load(file)
    
    # Extract tree from .jplace file
    tree = data['tree']

    # Extract placements
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            placements[gene_name] = placement['p'][0][1]  # Extract edge number
    
    return tree, placements

def find_closest_tip_labels(tree, placements):
    closest_tip_labels = {}
    for gene_id, edge_number in placements.items():
        closest_tip_labels[gene_id] = find_label_for_edge(tree, edge_number)
    return closest_tip_labels

def find_label_for_edge(tree, edge_number):
    # Implement the logic to find the tip label associated with the given edge number in the tree
    # You can use any method or library you prefer to parse the tree and find the label
    # Return "No matching label found" if no label is found for the given edge number
    pass

def print_closest_tip_labels(closest_tip_labels):
    for gene_id, tip_label in closest_tip_labels.items():
        print(f"Sample: {gene_id}, Closest Tip Label: {tip_label}")

def main(jplace_file):
    output_dir = './output'
    tree_path, _ = run_guppy(jplace_file, output_dir)
    tree, placements = extract_tree_and_placements(jplace_file)
    closest_tip_labels = find_closest_tip_labels(tree, placements)
    print_closest_tip_labels(closest_tip_labels)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python script.py <jplace_file>")
        sys.exit(1)
    main(sys.argv[1])
