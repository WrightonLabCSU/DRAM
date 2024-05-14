import re
from ete3 import Tree, TreeStyle, TextFace

def sanitize_newick(newick_file, sanitized_file):
    with open(newick_file, 'r') as f:
        newick_str = f.read()

    # Remove bootstrap values (e.g., [100])
    sanitized_str = re.sub(r'\[\d+\]', '', newick_str)

    with open(sanitized_file, 'w') as f:
        f.write(sanitized_str)

def color_labels(labels_file, newick_file, output_file, output_png):
    # Sanitize the Newick file
    sanitized_file = "sanitized_tree.nwk"
    sanitize_newick(newick_file, sanitized_file)

    # Read the labels from the text file
    with open(labels_file, 'r') as f:
        labels_to_color = set(line.strip() for line in f)

    # Parse the sanitized Newick file using ete3
    tree = Tree(sanitized_file)

    # Define the color you want to apply
    color = "red"

    # Traverse the tree and color the specified labels
    for leaf in tree.iter_leaves():
        if leaf.name in labels_to_color:
            face = TextFace(leaf.name, fgcolor=color)
            leaf.add_face(face, column=0, position="branch-right")

    # Set tree style to unrooted
    ts = TreeStyle()
    ts.mode = "c"
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 10  # Reduce the vertical space between branches
    ts.scale = 20  # Scale the tree

    # Render the tree to a PNG file
    tree.render(output_png, tree_style=ts)

    # Write the tree to a Newick file
    tree.write(outfile=output_file, format=1)

# File paths
labels_file = "labels.txt"
newick_file = "aligned_sequences.xml"
output_file = "colored_tree.nwk"
output_png = "colored_tree.png"

# Call the function
color_labels(labels_file, newick_file, output_file, output_png)

print('Colored labels and saved to colored_tree.png and colored_tree.nwk')
