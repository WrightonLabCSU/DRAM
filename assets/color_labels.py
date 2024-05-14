from Bio import Phylo
import sys

def color_labels(labels_file, newick_file, output_file, color='red'):
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

    print(f'Colored labels and saved to {output_file}')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file>")
        sys.exit(1)

    labels_file, newick_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]
    color_labels(labels_file, newick_file, output_file)
