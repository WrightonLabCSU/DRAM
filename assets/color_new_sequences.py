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

    # Function to recursively find all nodes and color the new sequences
    def color_nodes(element):
        for clade in element.findall('.//phy:clade', ns):
            name = clade.find('phy:name', ns)
            if name is not None and name.text in query_ids:
                color = ET.Element('phy:color', ns)
                red = ET.Element('phy:red', ns)
                green = ET.Element('phy:green', ns)
                blue = ET.Element('phy:blue', ns)
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
