import json
import pandas as pd
import sys
from Bio import Phylo

def load_phylogenetic_tree(tree_file):
    """Load the phylogenetic tree from a Newick file."""
    return Phylo.read(tree_file, 'newick')

def find_named_ancestor(tree, edge_number, tree_mapping):
    """Find the closest named ancestor in the tree from a given edge."""
    node = tree.find_clades({"name": str(edge_number)}).next()  # Assuming edge numbers are used as names
    if node is None:
        return None
    
    while node:
        if node.name in tree_mapping:
            return node.name, tree_mapping[node.name]
        node = node.parent
    
    return None

def extract_placements(jplace_path):
    """Extract placement edge numbers from the jplace file."""
    with open(jplace_path, 'r') as file:
        data = json.load(file)
    
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = placement['p'][0]
            edge_number = most_likely_placement[1]
            placements[gene_name] = edge_number
            print(f"Extracted placement for {gene_name}: Edge {edge_number}")
    return placements

def load_tree_mapping(mapping_path):
    mapping_df = pd.read_csv(mapping_path, sep='\t')
    tree_map = mapping_df.set_index('gene')['call'].to_dict()
    print(f"Loaded tree mapping for {len(tree_map)} genes with edges: {list(tree_map.keys())}")
    return tree_map

def update_annotations(annotations_path, placements, tree, tree_mapping):
    annotations_df = pd.read_csv(annotations_path, sep='\t')
    print(f"Attempting to update {len(annotations_df)} annotations.")
    
    for index, row in annotations_df.iterrows():
        gene_id = row['query_id']
        if gene_id in placements:
            edge = placements[gene_id]
            ancestor = find_named_ancestor(tree, edge, tree_mapping)
            if ancestor:
                annotations_df.at[index, 'tree-verified'] = ancestor[1]
                print(f"Gene {gene_id} placed on edge {edge} which corresponds to {ancestor[1]}")
            else:
                print(f"Gene {gene_id} placed on edge {edge} but no matching tree mapping found.")
        else:
            print(f"No placement found for gene {gene_id}")
    
    return annotations_df

def main(jplace_file, mapping_file, annotations_file, output_file):
    placements = extract_placements(jplace_file)
    tree_mapping = load_tree_mapping(mapping_file)
    updated_annotations = update_annotations(annotations_file, placements, tree_mapping)
    updated_annotations.to_csv(output_file, sep='\t', index=False)
    print("Updated annotations written to file.")
    
if __name__ == '__main__':
    jplace_file = sys.argv[1]
    mapping_file = sys.argv[2]
    annotations_file = sys.argv[3]
    output_file = sys.argv[4]
    tree_file = sys.argv[5]  # Add this line
    main(jplace_file, mapping_file, annotations_file, output_file, tree_file)
