import sys

def parse_prodigal_output(prod_out_file):
    gene_counter = {}  # Initialize a dictionary to keep track of gene counts per scaffold

    with open(prod_out_file, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'CDS':
                scaffold = parts[0]
                
                # Increment gene counter for the current scaffold
                if scaffold not in gene_counter:
                    gene_counter[scaffold] = 1
                else:
                    gene_counter[scaffold] += 1
                
                start = parts[3]
                stop = parts[4]
                
                # Append gene counter to scaffold name to create unique query_id
                query_id = f"{scaffold}_{gene_counter[scaffold]}"
                
                print(f"{query_id}\t{start}\t{stop}")

if __name__ == "__main__":
    prodigal_output_file = sys.argv[1]
    parse_prodigal_output(prodigal_output_file)
