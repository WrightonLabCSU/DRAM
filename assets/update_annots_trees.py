import json
import pandas as pd
import sys

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

def update_annotations(annotations_path, placements, tree_mapping):
    annotations_df = pd.read_csv(annotations_path, sep='\t')
    print(f"Attempting to update {len(annotations_df)} annotations.")
    
    missing_mappings = set()
    for index, row in annotations_df.iterrows():
        gene_id = row['query_id']
        if gene_id in placements:
            edge = placements[gene_id]
            if edge in tree_mapping:
                annotations_df.at[index, 'tree-verified'] = tree_mapping[edge]
                print(f"Gene {gene_id} placed on edge {edge} which corresponds to {tree_mapping[edge]}")
            else:
                missing_mappings.add(edge)
                print(f"Gene {gene_id} placed on edge {edge} but no matching tree mapping found.")
        else:
            print(f"No placement found for gene {gene_id}")

    print(f"Edges with missing mappings: {missing_mappings}")
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
    main(jplace_file, mapping_file, annotations_file, output_file)
