import json
import xml.etree.ElementTree as ET
import re
import sys

def parse_jplace(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data['tree'], jplace_data['placements']

def parse_tree(tree_newick):
    tree_newick = tree_newick.replace('{', '[').replace('}', ']')
    root = ET.Element("phyloxml", xmlns="http://www.phyloxml.org")
    phylogeny = ET.SubElement(root, "phylogeny", rooted="true")
    stack = [phylogeny]
    tokens = re.split(r"(\(|\)|,|;)", tree_newick)
    for token in tokens:
        token = token.strip()
        if not token:
            continue
        if token == '(':
            clade = ET.SubElement(stack[-1], "clade")
            stack.append(clade)
        elif token == ')':
            stack.pop()
        elif token == ',':
            stack.pop()
            clade = ET.SubElement(stack[-1], "clade")
            stack.append(clade)
        else:
            parts = token.split(":")
            name = parts[0]
            branch_length = parts[1] if len(parts) > 1 else "0.0"
            if name:
                ET.SubElement(stack[-1], "name").text = name
            if branch_length:
                ET.SubElement(stack[-1], "branch_length").text = branch_length.split('[')[0]
    return root

def insert_placements(root, placements):
    edge_map = defaultdict(list)
    for placement in placements:
        for p in placement['p']:
            edge_num = int(p[1])
            edge_map[edge_num].append(placement['nm'][0][0])

    def insert_into_tree(clade):
        for subclade in clade.findall("clade"):
            insert_into_tree(subclade)
        branch_length = clade.find("branch_length")
        if branch_length is not None:
            edge = re.search(r'\[(\d+)\]', branch_length.text)
            if edge and int(edge.group(1)) in edge_map:
                for seq_name in edge_map[int(edge.group(1))]:
                    new_clade = ET.SubElement(clade, "clade")
                    ET.SubElement(new_clade, "name").text = seq_name
                    ET.SubElement(new_clade, "branch_length").text = "0.0"

    insert_into_tree(root.find("phylogeny"))

def save_tree(root, output_path):
    tree = ET.ElementTree(root)
    ET.register_namespace("", "http://www.phyloxml.org")
    tree.write(output_path, encoding='utf-8', xml_declaration=True)

def main():
    if len(sys.argv) != 4:
        print("Usage: python update_tree.py <jplace_path> <input_xml_path> <output_xml_path>")
        sys.exit(1)

    jplace_path, input_xml_path, output_xml_path = sys.argv[1:]

    # Parse the jplace file
    tree_newick, placements = parse_jplace(jplace_path)

    # Parse the input XML tree
    tree_root = parse_tree(tree_newick)

    # Insert placements into the tree
    insert_placements(tree_root, placements)

    # Save the updated tree to an output XML file
    save_tree(tree_root, output_xml_path)

if __name__ == "__main__":
    main()
