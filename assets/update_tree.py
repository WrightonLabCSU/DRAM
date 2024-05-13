import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

def parse_jplace(jplace_path):
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data['tree'], jplace_data['placements']

def parse_tree(tree_newick):
    tree_newick = tree_newick.replace('{', '[').replace('}', ']')
    root = ET.Element("phy:phylogeny", rooted="true")
    tree_element = ET.SubElement(root, "phy:clade")
    stack = [tree_element]
    for token in re.split(r"(\(|\)|,|;)", tree_newick):
        token = token.strip()
        if not token:
            continue
        if token == '(':
            clade = ET.SubElement(stack[-1], "phy:clade")
            stack.append(clade)
        elif token == ')':
            stack.pop()
        elif token == ',':
            stack.pop()
            clade = ET.SubElement(stack[-1], "phy:clade")
            stack.append(clade)
        else:
            name, branch_length = token.split(":")
            stack[-1].append(ET.Element("phy:name", text=name))
            stack[-1].append(ET.Element("phy:branch_length", text=branch_length))
    return root

def insert_placements(root, placements):
    edge_map = defaultdict(list)
    for placement in placements:
        for p in placement['p']:
            edge_num = p[1]
            edge_map[edge_num].append(placement['nm'][0][0])

    def insert_into_tree(clade):
        for subclade in clade.findall("phy:clade", namespaces=ns):
            insert_into_tree(subclade)
        edge = int(clade.find("phy:branch_length", namespaces=ns).text.split('{')[1].split('}')[0])
        if edge in edge_map:
            for seq_name in edge_map[edge]:
                new_clade = ET.SubElement(clade, "phy:clade")
                ET.SubElement(new_clade, "phy:name").text = seq_name
                ET.SubElement(new_clade, "phy:branch_length").text = "0.0"

    insert_into_tree(root.find("phy:clade", namespaces=ns))

def save_tree(root, output_path):
    tree = ET.ElementTree(root)
    ET.register_namespace("phy", "http://www.phyloxml.org")
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
