import pandas as pd
import sys

def load_tree_mapping(mapping_tsv):
    df = pd.read_csv(mapping_tsv, sep='\t')
    return dict(zip(df['gene'], df['call']))

def parse_placement_csv(csv_path, tree_mapping):
    df = pd.read_csv(csv_path)
    placement_map = {}
    for _, row in df.iterrows():
        read_name = row['name']
        closest_leaf = row['edge_num']
        mapped_value = tree_mapping.get(closest_leaf, "No mapping found")
        placement_map[read_name] = f"{mapped_value};{closest_leaf}"
        print(f"{read_name} classified as {mapped_value} (Closest Leaf: {closest_leaf})")
    return placement_map

def update_tsv(tsv_path, output_tsv_path, placement_map):
    df = pd.read_csv(tsv_path, sep='\t')
    df['tree_verified'] = df['query_id'].map(placement_map).fillna('')

    # Reorder columns to place 'tree_verified' after 'gene_number'
    col_order_start = df.columns.tolist()[:df.columns.get_loc('gene_number')+1] + ['tree_verified']
    col_order_end = [col for col in df.columns if col not in col_order_start]
    df = df[col_order_start + col_order_end]

    df.to_csv(output_tsv_path, sep='\t', index=False)

def main():
    if len(sys.argv) != 5:
        print("Usage: python update_annots_trees.py <tsv_path> <mapping_tsv> <csv_path> <output_tsv_path>")
        sys.exit(1)

    tsv_path, mapping_tsv, csv_path, output_tsv_path = sys.argv[1:]
    tree_mapping = load_tree_mapping(mapping_tsv)
    placement_map = parse_placement_csv(csv_path, tree_mapping)
    update_tsv(tsv_path, output_tsv_path, placement_map)

if __name__ == "__main__":
    main()
