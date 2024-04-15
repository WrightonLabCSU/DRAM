import json
import sys

def load_jplace_file(jplace_path):
    """ Load and parse the .jplace JSON file. """
    with open(jplace_path, 'r') as file:
        jplace_data = json.load(file)
    return jplace_data

def extract_placement_details(jplace_data):
    """ Extract and print details from placements. """
    placements = jplace_data['placements']
    results = []
    for placement in placements:
        for placement_detail in placement['p']:
            edge_num = placement_detail[1]
            likelihood = placement_detail[3]
            for name, _ in placement['nm']:
                results.append((name, edge_num, likelihood))
    return results

def print_placements(results):
    """ Print the placement details for each sequence. """
    for name, edge_num, likelihood in results:
        print(f"Sequence: {name}, Edge: {edge_num}, Likelihood: {likelihood}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python update_annots_trees.py <jplace_path>")
        return
    
    jplace_path = sys.argv[1]
    jplace_data = load_jplace_file(jplace_path)
    results = extract_placement_details(jplace_data)
    print_placements(results)

if __name__ == "__main__":
    main()
