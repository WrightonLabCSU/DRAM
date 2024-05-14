import sys
from ete3 import Tree, NodeStyle, TreeStyle, AttrFace, faces

def color_labels(labels_file, newick_file, output_file, output_png):
    # Read the labels from the text file
    with open(labels_file, 'r') as f:
        labels_to_color = set(line.strip() for line in f)

    # Load the tree
    tree = Tree(newick_file)

    # Color the labels
    for leaf in tree:
        if leaf.name in labels_to_color:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "red"
            leaf.set_style(nstyle)

    # Save the modified tree to a new Newick file
    tree.write(outfile=output_file, format=1)

    # Create and configure tree style
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True

    # Render the tree to a PNG file
    tree.render(output_png, tree_style=ts)

    print('Colored labels and saved to', output_png)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file> <output_png>")
        sys.exit(1)
    
    labels_file = sys.argv[1]
    newick_file = sys.argv[2]
    output_file = sys.argv[3]
    output_png = sys.argv[4]
    
    color_labels(labels_file, newick_file, output_file, output_png)
