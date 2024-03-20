import argparse
import csv
from collections import defaultdict
import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate GFF and/or GBK files from raw annotations, with specified databases formatting.")
    parser.add_argument("--gff", action='store_true', help="Generate GFF file")
    parser.add_argument("--gbk", action='store_true', help="Generate GBK file")
    parser.add_argument("--samples_paths", nargs='+', help="Alternating list of sample names and paths to their .fna files.")
    parser.add_argument("--database_list", type=str, help="List of databases to include in the annotations. Use 'empty' for all.", default="empty")
    parser.add_argument("--annotations", required=True, help="Path to the raw annotations file")
    args = parser.parse_args()
    args.database_list = None if args.database_list == "empty" else args.database_list.split()
    return args

import argparse
import csv
from collections import defaultdict
import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
import logging
import urllib.parse  # For escaping characters in accordance with RFC 3986

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate GFF and/or GBK files from raw annotations, with specified databases formatting.")
    parser.add_argument("--gff", action='store_true', help="Generate GFF file")
    parser.add_argument("--gbk", action='store_true', help="Generate GBK file")
    parser.add_argument("--samples_paths", nargs='+', help="Alternating list of sample names and paths to their .fna files.")
    parser.add_argument("--database_list", type=str, help="Comma-separated list of databases to include in the annotations. Use 'empty' for all.", default="empty")
    parser.add_argument("--annotations", required=True, help="Path to the raw annotations file")
    args = parser.parse_args()
    args.database_list = None if args.database_list == "empty" else args.database_list.split(',')
    return args

def escape_gff3_value(value):
    """
    Escapes characters in a string according to GFF3 specification using URL encoding.
    """
    return urllib.parse.quote(value, safe=':/,=.-')

def parse_samples_and_paths(samples_paths):
    """
    Parses the provided list of sample names and .fna file paths into a structured dictionary.
    """
    cleaned_samples_paths = [item.strip("[]',") for item in samples_paths]
    iterator = iter(cleaned_samples_paths)
    return dict(zip(iterator, iterator))

def format_attributes(annotation, database_list):
    """
    Format and order database-specific annotations for the GFF attributes column, with customized formatting and escaping.
    """
    attributes = []
    for key, value in sorted(annotation.items()):
        if key.endswith('_id') or key.endswith('_description'):
            db_name = key.split('_')[0]
            if database_list is None or db_name in database_list:
                upper_db_name = db_name.upper()
                if key.endswith('_id'):
                    description_key = f"{db_name}_description"
                    desc = escape_gff3_value(annotation.get(description_key, "NA"))
                    attributes.append(f"{upper_db_name}_ID={value};{upper_db_name}_Description={desc}")
    return "; ".join(attributes)

def generate_gff(samples_annotations, database_list):
    """Generate GFF files for each sample, filtered by specified databases."""
    os.makedirs("GFF", exist_ok=True)
    for sample, annotations in samples_annotations.items():
        with open(f"GFF/{sample}.gff", "w") as gff_file:
            gff_file.write("##gff-version 3\n")
            for annotation in annotations:
                # Ensure seqid complies with GFF3 specification by escaping reserved characters
                seqid = escape_gff3_value(annotation['query_id'])
                source = '.'  # Placeholder, should be replaced with actual source if available
                type = "gene"  # Example type, should be replaced with actual type from annotation if available
                start = annotation['start_position']
                end = annotation['stop_position']
                score = '.'  # Placeholder, replace if score available
                strand = '+' if annotation['strandedness'] == '+1' else '-'
                phase = '.'  # Placeholder, applicable for CDS features
                attributes_str = format_attributes(annotation, database_list)
                gff_line = f"{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}\n"
                gff_file.write(gff_line)

def main():
    args = parse_arguments()

    # Load annotations and organize by sample
    samples_annotations = defaultdict(list)
    with open(args.annotations, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            samples_annotations[row['sample']].append(row)

    # Directly parse the samples and paths passed as arguments
    samples_and_paths = parse_samples_and_paths(args.samples_paths)

def parse_fna_sequence(fna_file_path):
    """Parse the .fna file to get sequences indexed by their header name."""
    sequences = {}
    for seq_record in SeqIO.parse(fna_file_path, "fasta"):
        # The header name is directly used as the key
        header_name = seq_record.id
        sequences[header_name] = seq_record.seq
    return sequences

def format_qualifiers(annotation, database_list):
    """
    Format database-specific annotations into qualifiers.
    """
    qualifiers = {}
    for key, value in annotation.items():
        if key.endswith('_id') or key.endswith('_description'):
            db_name = key.split('_')[0]
            upper_db_name = db_name.upper()
            if database_list is None or db_name in database_list:
                description_key = f"{db_name}_description"
                desc = annotation.get(description_key, "NA")
                qualifiers[f"({upper_db_name}) {db_name}_id"] = value
                qualifiers[f"({upper_db_name}) description"] = desc
    return qualifiers

def find_sample_fna_files(sample, fna_directory):
    """
    Find all .fna files starting with the sample name in the specified directory.
    """
    pattern = os.path.join(fna_directory, f"{sample}*.fna")
    return glob.glob(pattern)

def aggregate_sample_sequences(sample_files):
    """
    Aggregate sequences from multiple .fna files for a sample.
    """
    sequences = {}
    for file_path in sample_files:
        for seq_record in SeqIO.parse(file_path, "fasta"):
            sequences[seq_record.id] = seq_record.seq
    return sequences

def generate_gbk(samples_annotations, database_list, samples_and_paths):
    print("Starting GBK generation...")
    os.makedirs("GBK", exist_ok=True)

    print(f"Cleaned and parsed samples and paths: {samples_and_paths}")

    for sample, annotations in samples_annotations.items():
        print(f"\nProcessing sample: {sample}")

        if sample in samples_and_paths:
            fna_file_path = samples_and_paths[sample]
            print(f"Using .fna file for {sample}: {fna_file_path}")

            if os.path.exists(fna_file_path):
                sequences = parse_fna_sequence(fna_file_path)

                # Initialize an empty SeqRecord with basic info
                seq_record = SeqRecord(Seq(""), id=sample, description=f"Generated GBK file for {sample}", annotations={"molecule_type": "DNA", "source": "DRAM2"})
                
                if annotations:
                    metadata = annotations[0]  # Assuming shared metadata across each sample's annotations
                    taxonomy_info = metadata.get('taxonomy', 'Not Available')
                    seq_record.annotations["organism"] = taxonomy_info
                    completeness_info = str(metadata.get('Completeness', 'Not Available'))
                    contamination_info = str(metadata.get('Contamination', 'Not Available'))
                    seq_record.annotations["note"] = f"Completeness: {completeness_info}; Contamination: {contamination_info}"

                for annotation in annotations:
                    query_id = annotation['query_id']
                    if query_id in sequences:
                        # Retrieve sequence for this annotation
                        sequence = sequences[query_id]
                        feature_location = FeatureLocation(start=int(annotation['start_position']) - 1, end=int(annotation['stop_position']), strand=1 if annotation['strandedness'] == '+1' else -1)
                        qualifiers = format_qualifiers(annotation, database_list)
                        feature = SeqFeature(feature_location, type="gene", qualifiers=qualifiers)
                        seq_record.features.append(feature)

                        # Concatenate sequence to SeqRecord.seq for the entire sample
                        # This assumes sequences are non-overlapping and can be concatenated
                        seq_record.seq += sequence

                output_filename = f"GBK/{sample}.gbk"
                with open(output_filename, "w") as output_handle:
                    SeqIO.write([seq_record], output_handle, "genbank")
                print(f"GBK file generated for {sample}: {output_filename}")
            else:
                print(f"File does not exist: {fna_file_path}")
        else:
            print(f"No .fna file path found for sample {sample}")


def main():
    args = parse_arguments()

    # Load annotations and organize by sample
    samples_annotations = defaultdict(list)
    with open(args.annotations, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            samples_annotations[row['sample']].append(row)

    # Directly parse the samples and paths passed as arguments
    samples_and_paths = parse_samples_and_paths(args.samples_paths)

    if args.gbk:
        generate_gbk(samples_annotations, args.database_list, samples_and_paths)

    if args.gff:
        generate_gff(samples_annotations, args.database_list)


if __name__ == "__main__":
    main()
