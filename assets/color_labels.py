from Bio import Phylo
import sys
import matplotlib.pyplot as plt
import networkx as nx

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
        if clade.name is None:
            clade.name = f"Unnamed_{id(clade)}"

    # Write the modified tree to a new file
    Phylo.write(tree, output_file, 'newick')

    # Convert tree to a NetworkX graph
    def to_networkx(tree, attrs={}):
        """Convert a Biopython tree to a NetworkX graph."""
        graph = nx.DiGraph()
        def add_edge(parent, child):
            graph.add_edge(parent.name, child.name, **attrs)
            for clade in child.clades:
                add_edge(child, clade)
        for clade in tree.clade.clades:
            add_edge(tree.clade, clade)
        return graph

    # Draw the unrooted tree using NetworkX
    graph = to_networkx(tree)
    pos = nx.spring_layout(graph, k=0.5, iterations=50)
    plt.figure(figsize=(20, 20))  # Adjust the figure size as needed
    nx.draw(graph, pos, with_labels=True, font_weight='bold', font_size=10, node_size=0, node_color='w', edge_color='k')

    # Save the plot as a PNG file
    plt.savefig(output_png)
    plt.close()

    print(f'Colored labels saved to {output_file}')
    print(f'Unrooted tree saved as PNG to {output_png}')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python color_labels.py <labels_file> <newick_file> <output_file> <output_png>")
        sys.exit(1)

    labels_file, newick_file, output_file, output_png = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    color_labels(labels_file, newick_file, output_file, output_png)
