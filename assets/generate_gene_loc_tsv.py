import sys

def parse_prodigal_output(prod_out_file):
    with open(prod_out_file, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'CDS':
                attributes = parts[8].split(';')
                id_attr = [attr for attr in attributes if attr.startswith('ID=')][0]
                query_id = id_attr.split('=')[1]  # Extract the ID value
                
                start = parts[3]
                stop = parts[4]
                
                print(f"{query_id}\t{start}\t{stop}")

if __name__ == "__main__":
    prodigal_output_file = sys.argv[1]
    parse_prodigal_output(prodigal_output_file)
