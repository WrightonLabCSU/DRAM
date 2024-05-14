import xml.etree.ElementTree as ET
import sys

def color_labels(labels_file, xml_file, output_file, color='red'):
    # Read the labels from the text file
    with open(labels_file, 'r') as f:
        labels_to_color = set(line.strip() for line in f)

    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Function to add color attribute to a label
    def color_label(label):
        label.set('color', color)

    # Traverse the XML tree and color the specified labels
    for elem in root.iter():
        if elem.tag == 'name' and elem.text in labels_to_color:
            parent = elem.getparent()
            for label in parent.iter('label'):
                color_label(label)

    # Save the modified XML to a new file
    tree.write(output_file)

    print(f'Colored labels and saved to {output_file}')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python color_labels.py <labels_file> <xml_file> <output_file>")
        sys.exit(1)

    labels_file, xml_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    color_labels(labels_file, xml_file, output_file)
