import json
import xml.etree.ElementTree as ET
import re
from collections import defaultdict
import sys
from ete3 import Tree, Phyloxml

def parse_jplace(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data['tree'], jplace_data['placements']

def add_placements_to_tree(tree_newick, placements):
    # Load the tree
    tree = Tree(tree_newick, format=1)
    
    # Add placements
    for placement in placements:
        for p in placement['p']:
            edge_num = p[1]
            name = placement['nm'][0][0]
            edge = tree.search_nodes(name=str(edge_num))
            if edge:
                new_clade = Tree(name=name, dist=0.0)
                edge[0].add_child(new_clade)
    
    return tree

def save_tree_as_phyloxml(tree, output_path):
    phyloxml = tree.write(format=5)
    with open(output_path, 'w') as file:
        file.write(phyloxml)

def main():
    jplace_path = 'aligned_sequences.jplace'  # Replace with the actual path if different
    output_xml_path = 'aligned_sequences_updated.xml'  # Replace with the actual path if different

    # Parse the jplace file
    tree_newick, placements = parse_jplace(jplace_path)

    # Add placements to the tree
    tree = add_placements_to_tree(tree_newick, placements)

    # Save the tree in phyloXML format
    save_tree_as_phyloxml(tree, output_xml_path)

if __name__ == "__main__":
    main()
