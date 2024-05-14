import sys
from ete3 import Tree, TreeStyle, NodeStyle

def color_labels(labels_file, newick_file, output_file, output_png):
    # Read the labels from the text file
    with open(labels_file, 'r') as f:
        labels_to_color = set(line.strip() for line in f)

    # Load the tree from the newick file
    tree = Tree(newick_file)

    # Define the color you want to apply
    color = "red"

    # Traverse the tree and color the specified labels
    for leaf in tree.iter_leaves():
        if leaf.name in labels_to_color:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color
            leaf.set_style(nstyle)
    
    # Save the modified tree to a new newick file
    tree.write(outfile=output_file)

    # Visualize the tree and save as PNG
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = True
    ts.scale = 20
    ts.branch_vertical_margin = 10

    tree.render(output_png, tree_style=ts)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file> <output_png>")
        sys.exit(1)

    labels_file, newick_file, output_file, output_png = sys.argv[1:]
    color_labels(labels_file, newick_file, output_file, output_png)
