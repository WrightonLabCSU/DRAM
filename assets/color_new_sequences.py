import sys
import xml.etree.ElementTree as ET

def color_new_sequences(xml_path, query_ids_path, output_path):
    # Load the list of query IDs to be colored
    with open(query_ids_path, 'r') as f:
        query_ids = set(line.strip().split('\t')[1] for line in f.readlines())

    # Parse the XML tree
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # Function to recursively find all nodes and color the new sequences
    def color_nodes(element):
        for clade in element.findall('.//{http://www.phyloxml.org}clade'):
            name = clade.find('{http://www.phyloxml.org}name')
            if name is not None and name.text in query_ids:
                color = ET.Element('{http://www.phyloxml.org}color')
                red = ET.Element('{http://www.phyloxml.org}red')
                green = ET.Element('{http://www.phyloxml.org}green')
                blue = ET.Element('{http://www.phyloxml.org}blue')
                red.text = '255'
                green.text = '0'
                blue.text = '0'
                color.append(red)
                color.append(green)
                color.append(blue)
                clade.append(color)
            # Recursively color nested clades
            color_nodes(clade)

    # Start coloring from the root element
    color_nodes(root)

    # Write the modified tree to the output file
    tree.write(output_path)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python color_new_sequences.py <input_xml> <query_ids> <output_xml>")
        sys.exit(1)

    input_xml, query_ids, output_xml = sys.argv[1:]
    color_new_sequences(input_xml, query_ids, output_xml)
