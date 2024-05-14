import re
import sys
from ete3 import Tree, TreeStyle, NodeStyle

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
    tree = Tree(cleaned_newick_file)

    # Define the color you want to apply
    color = "red"

    # Traverse the tree and color the specified labels
    for leaf in tree.iter_leaves():
        if leaf.name in labels_to_color:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color
            leaf.set_style(nstyle)
    
    # Save the modified tree to a new Newick file
    tree.write(outfile=output_file)

    # Visualize the tree and save as SVG and PNG
    ts = TreeStyle()
    ts.mode = "c"  # Circular mode
    ts.show_leaf_name = True
    ts.scale = 200  # Increase scale for larger image
    ts.branch_vertical_margin = 10

    # Render the tree to SVG and PNG files without opening a display
    tree.render(output_svg, tree_style=ts)
    tree.render(output_png, tree_style=ts, dpi=300)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file> <output_svg> <output_png>")
        sys.exit(1)

    labels_file, newick_file, output_file, output_svg, output_png = sys.argv[1:]
    color_labels(labels_file, newick_file, output_file, output_svg, output_png)
