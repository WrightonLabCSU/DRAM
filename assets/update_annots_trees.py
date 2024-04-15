import json
import subprocess
import os

def run_guppy(jplace_file, output_dir):
    subprocess.run(['guppy', 'tog', jplace_file, '-o', f"{output_dir}/tree_with_placements.newick"], check=True)
    subprocess.run(['guppy', 'edpl', '--csv', jplace_file, '-o', f"{output_dir}/edpl.csv"], check=True)
    return f"{output_dir}/tree_with_placements.newick", f"{output_dir}/edpl.csv"

def extract_placements(jplace_file):
    with open(jplace_file, 'r') as file:
        data = json.load(file)
    placements = {}
    for placement in data['placements']:
        for gene_info in placement['nm']:
            gene_name, _ = gene_info
            most_likely_placement = min(placement['p'], key=lambda x: x[3])  # Minimize by likelihood score
            edge_number, _, _ = most_likely_placement[1], most_likely_placement[2], most_likely_placement[3]
            placements[gene_name] = edge_number
    return placements

def print_closest_tip_labels(placements):
    for gene_id, edge_number in placements.items():
        print(f"Sample: {gene_id}, Closest Edge Number: {edge_number}")

def main(jplace_file):
    output_dir = './output'
    os.makedirs(output_dir, exist_ok=True)
    _, _ = run_guppy(jplace_file, output_dir)
    placements = extract_placements(jplace_file)
    print_closest_tip_labels(placements)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script.py <jplace_file>")
        sys.exit(1)
    main(sys.argv[1])
