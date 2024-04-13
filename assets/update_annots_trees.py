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
    return placements

def load_tree_mapping(mapping_path):
    """Load tree mapping from a TSV file."""
    mapping_df = pd.read_csv(mapping_path, sep='\t')
    return mapping_df.set_index('gene')['call'].to_dict()

def update_annotations(annotations_path, placements, tree_mapping):
    """Update the annotations with tree-verified subunit information."""
    annotations_df = pd.read_csv(annotations_path, sep='\t')
    annotations_df['tree-verified'] = annotations_df['query_id'].map(placements).map(tree_mapping)
    return annotations_df

def main(jplace_file, mapping_file, annotations_file, output_file):
    placements = extract_placements(jplace_file)
    tree_mapping = load_tree_mapping(mapping_file)
    updated_annotations = update_annotations(annotations_file, placements, tree_mapping)
    updated_annotations.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    jplace_file = sys.argv[1]
    mapping_file = sys.argv[2]
    annotations_file = sys.argv[3]
    output_file = sys.argv[4]
    main(jplace_file, mapping_file, annotations_file, output_file)
