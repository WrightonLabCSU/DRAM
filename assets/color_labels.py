import re
import sys
import matplotlib.pyplot as plt
from Bio import Phylo
from ete3 import Tree

def clean_newick(newick_file, cleaned_newick_file):
    with open(newick_file, 'r') as infile:
        newick_data = infile.read()
    cleaned_data = re.sub(r'\[\d+\]', '', newick_data)
    with open(cleaned_newick_file, 'w') as outfile:
        outfile.write(cleaned_data)

def color_labels(labels_file, newick_file, output_file, output_svg, output_png):
    # Clean the Newick file
    cleaned_newick_file = "cleaned_" + newick_file
    clean_newick(newick_file, cleaned_newick_file)
    
    # Read the labels from the text file
    with open(labels_file, 'r') as f:
        labels_to_color = set(line.strip() for line in f)

    # Load the cleaned tree from the Newick file
    tree = Phylo.read(cleaned_newick_file, "newick")

    # Define the color you want to apply
    color = "red"

    # Traverse the tree and color the specified labels
    for clade in tree.find_clades():
        if clade.name in labels_to_color:
            clade.color = color
    
    # Save the modified tree to a new Newick file
    Phylo.write(tree, output_file, "newick")

    # Visualize the tree and save as SVG and PNG
    fig = plt.figure(figsize=(15, 15))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig(output_svg, format="svg")
    plt.savefig(output_png, format="png", dpi=300)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file> <output_svg> <output_png>")
        sys.exit(1)

    labels_file, newick_file, output_file, output_svg, output_png = sys.argv[1:]
    color_labels(labels_file, newick_file, output_file, output_svg, output_png)
