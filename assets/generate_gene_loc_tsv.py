# generate_gene_loc_tsv.py
import sys

def parse_prodigal_output(prod_out_file):
    with open(prod_out_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'CDS':
                scaffold = parts[0]
                start = parts[3]
                stop = parts[4]
                print(f"{scaffold}\t{start}\t{stop}")

if __name__ == "__main__":
    prodigal_output_file = sys.argv[1]
    parse_prodigal_output(prodigal_output_file)
