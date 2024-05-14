from Bio import Phylo
import sys
import matplotlib.pyplot as plt

def color_labels(labels_file, newick_file, output_file, output_png, color='red'):
    # Read the labels from the text file
    with open(labels_file, 'r') as f:
        labels_to_color = set(line.strip() for line in f)

    # Parse the Newick tree
    tree = Phylo.read(newick_file, 'newick')

    # Traverse the tree and color the specified labels
    for clade in tree.find_clades():
        if clade.name in labels_to_color:
            clade.name = f"[&!color={color}]{clade.name}"

    # Write the modified tree to a new file
    Phylo.write(tree, output_file, 'newick')

    # Draw the tree as an unrooted tree and save as PNG
    fig = plt.figure(figsize=(10, 10))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=axes, branch_labels=None)
    fig.savefig(output_png)
    plt.close(fig)

    print(f'Colored labels saved to {output_file}')
    print(f'Unrooted tree saved as PNG to {output_png}')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file> <output_png>")
        sys.exit(1)

    labels_file, newick_file, output_file, output_png = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    color_labels(labels_file, newick_file, output_file, output_png)
