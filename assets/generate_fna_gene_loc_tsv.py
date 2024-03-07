import sys

def parse_fna_file(input_file):
    """
    Parse the input .fna file and extract the query name, start position,
    and stop position for each entry.

    Parameters:
    input_file (str): The path to the .fna file to be parsed.

    Returns:
    list of tuples: Each tuple contains (query_name, start_position, stop_position).
    """
    gene_locations = []
    with open(input_file, 'r') as fna:
        for line in fna:
            if line.startswith('>'):
                parts = line.split('#')
                if len(parts) >= 4:
                    query_name = parts[0].strip('>')
                    start_pos = parts[1].strip()
                    stop_pos = parts[2].strip()
                    gene_locations.append((query_name, start_pos, stop_pos))
    return gene_locations

def write_tsv(output_file, gene_locations):
    """
    Write the extracted gene locations to a TSV file.

    Parameters:
    output_file (str): The path to the output .tsv file.
    gene_locations (list of tuples): Each tuple contains (query_name, start_position, stop_position).
    """
    with open(output_file, 'w') as tsv:
        for query_name, start_pos, stop_pos in gene_locations:
            tsv.write(f"{query_name}\t{start_pos}\t{stop_pos}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_fna_gene_loc_tsv.py <input_fna_file> <output_tsv_file>")
        sys.exit(1)
    
    input_fna_file = sys.argv[1]
    output_tsv_file = sys.argv[2]
    
    # Parse the .fna file to get gene locations
    gene_locations = parse_fna_file(input_fna_file)
    
    # Write the gene locations to a .tsv file
    write_tsv(output_tsv_file, gene_locations)
