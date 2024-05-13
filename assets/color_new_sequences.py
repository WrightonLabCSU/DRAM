import sys
import xml.etree.ElementTree as ET

def color_new_sequences(xml_path, query_ids_path, output_path):
    # Load the list of query IDs to be colored
    with open(query_ids_path, 'r') as f:
        query_ids = set(line.strip().split('\t')[1] for line in f.readlines())

    # Parse the XML tree
    tree = ET.parse(xml_path)
    root = tree.getroot()

    # Define the phyloXML namespace
    ns = {'phy': 'http://www.phyloxml.org'}

    # Function to color the labels of new sequences
    def color_labels(element):
        for clade in element.findall('.//phy:clade', ns):
            name = clade.find('phy:name', ns)
            if name is not None and name.text in query_ids:
                # Create a new color element for the sequence label
                color_elem = ET.Element("phy:color", {'r': '255', 'g': '0', 'b': '0'})
                name.append(color_elem)
            # Recursively color nested clades
            color_labels(clade)

    # Start coloring from the root element
    color_labels(root)

    # Register namespace
    ET.register_namespace('', "http://www.phyloxml.org")

    # Write the modified tree to the output file
    tree.write(output_path, encoding='utf-8', xml_declaration=True)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python color_new_sequences.py <input_xml> <query_ids> <output_xml>")
        sys.exit(1)

    input_xml, query_ids, output_xml = sys.argv[1:]
    color_new_sequences(input_xml, query_ids, output_xml)
