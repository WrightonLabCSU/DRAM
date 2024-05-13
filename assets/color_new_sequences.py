import xml.etree.ElementTree as ET
import sys

def color_new_sequences(xml_path, new_sequences, output_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    for node in root.findall('.//clade'):
        name_elem = node.find('name')
        if name_elem is not None and name_elem.text in new_sequences:
            branch_color = ET.SubElement(node, 'branch_color')
            branch_color.set('r', '255')
            branch_color.set('g', '0')
            branch_color.set('b', '0')
            branch_color.set('a', '1.0')

    tree.write(output_path)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python color_new_sequences.py <xml_path> <new_sequences_file> <output_path>")
        sys.exit(1)

    xml_path = sys.argv[1]
    new_sequences_file = sys.argv[2]
    output_path = sys.argv[3]

    with open(new_sequences_file) as f:
        new_sequences = [line.strip() for line in f]

    color_new_sequences(xml_path, new_sequences, output_path)
