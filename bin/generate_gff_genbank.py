import argparse
import csv
from collections import defaultdict
import os
import glob
from datetime import datetime
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

def parse_samples_and_paths(samples_paths):
    """
    Parses the provided list of sample names and .fna file paths into a structured dictionary.
    """
    cleaned_samples_paths = [item.strip("[]',") for item in samples_paths]
    iterator = iter(cleaned_samples_paths)
    return dict(zip(iterator, iterator))

def sanitize_taxonomy(taxonomy):
    """
    Sanitizes the taxonomy string for safe inclusion in GFF comments.
    Replaces semicolons with commas or another safe delimiter to avoid parsing issues.
    """
    return taxonomy.replace(';', ',')  # Replace semicolons with commas

def escape_gff3_value(value):
    """
    Escapes characters in a string for GFF3 compliance, using URL encoding.
    """
    # RFC 3986 percent-encoding, excluding reserved characters used in GFF3 format
    return urllib.parse.quote(value, safe=':/,=.-')

def format_attributes(annotation, database_list):
    """
    Format and order database-specific annotations for the GFF attributes column, with customized formatting.
    This version omits empty attributes and ensures no extraneous spaces or semicolons.
    """
    attributes = []
    for key, value in sorted(annotation.items()):
        if (key.endswith('_id') or key.endswith('_description')) and value:
            db_name = key.split('_')[0]
            if database_list is None or db_name in database_list:
                # Ensure value is URL-encoded
                encoded_value = escape_gff3_value(value)
                attributes.append(f"{db_name.upper()}_{key.upper()}={encoded_value}")
    
    # Join attributes with semicolons, avoiding spaces for GFF3 compliance
    return ";".join(attributes)

def generate_gff(samples_annotations, database_list):
    """
    Generate GFF files for each sample, filtered by specified databases,
    ensuring adherence to GFF3 specifications, including attribute encoding and metadata comments.
    """
    os.makedirs("GFF", exist_ok=True)
    for sample, annotations in samples_annotations.items():
        with open(f"GFF/{sample}.gff", "w") as gff_file:
            gff_file.write("##gff-version 3\n")
            
            # Extract metadata from the first annotation (assuming it's consistent across annotations)
            metadata = annotations[0]
            completeness = metadata.get('Completeness', 'NA')
            contamination = metadata.get('Contamination', 'NA')
            taxonomy = escape_gff3_value(metadata.get('taxonomy', 'NA').replace(';', ','))
            
            gff_file.write(f"# Completeness: {completeness}\n")
            gff_file.write(f"# Contamination: {contamination}\n")
            gff_file.write(f"# Taxonomy: {taxonomy}\n")
            
            for annotation in annotations:
                try:
                    seqid = escape_gff3_value(annotation['query_id'])
                    source = '.'
                    type = "gene"
                    start = annotation['start_position']
                    end = annotation['stop_position']
                    score = '.'
                    strand = '+' if annotation['strandedness'] == '+1' else '-'
                    phase = '.'
                    attributes_str = format_attributes(annotation, database_list)
                    
                    gff_line = f"{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}\n"
                    gff_file.write(gff_line)
                except KeyError as e:
                    logging.warning(f"Missing key in annotation: {e}. Skipping annotation: {annotation}")



def parse_fna_sequence(fna_file_path):
    """Parse the .fna file to get sequences indexed by their header name."""
    sequences = {}
    for seq_record in SeqIO.parse(fna_file_path, "fasta"):
        header_name = seq_record.id
        sequences[header_name] = seq_record.seq
    return sequences

def aggregate_sample_sequences(samples_and_paths):
    """
    Aggregate sequences from .fna files, ensuring correct matching with sample names.
    This function adjusts to the fact that each .fna file might contain multiple sequences.
    """
    sequences = {}
    for sample, file_path in samples_and_paths.items():
        if os.path.exists(file_path):
            # Assuming we're collecting all sequences for a sample in a list
            sequences[sample] = [seq_record.seq for seq_record in SeqIO.parse(file_path, "fasta")]
        else:
            print(f"File does not exist: {file_path}")
    return sequences

def format_qualifiers(annotation, database_list=None):
    """
    Formats qualifiers for a GenBank feature, ensuring correct syntax and inclusion of essential information.
    Optionally filters the qualifiers based on a list of specified databases.
    """
    qualifiers = {}
    for key, value in annotation.items():
        if key in ['gene', 'product', 'note', 'function', 'locus_tag']:  # Add 'locus_tag' if relevant
            qualifiers[key] = value
        else:
            # Check if the key is database-specific and whether it should be included
            db_name = key.split('_')[0]
            if database_list is None or db_name in database_list:
                qualifiers[key] = value
    return qualifiers

def generate_gbk_feature(annotation):
    """Generate a SeqFeature for a GenBank file based on annotation data."""
    start = int(annotation['start_position']) - 1  # Convert to 0-based indexing
    end = int(annotation['stop_position'])  # End is inclusive, no adjustment needed
    strand = 1 if annotation['strandedness'] == '+1' else -1  # Convert strandedness
    location = FeatureLocation(start, end, strand=strand)

    qualifiers = {
        'locus_tag': annotation.get('query_id', 'unknown_locus'),
        # Assuming 'product' or a general description can be derived from the data
        'product': next((annotation.get(k) for k in annotation.keys() if k.endswith('_description') and annotation[k]), 'unknown product')
    }

    # Dynamically add database-specific information as qualifiers
    for key, value in annotation.items():
        if key.endswith(('_id', '_description', '_EC')) and value:
            # Simple heuristic to extract database prefix and adjust it if necessary
            db_prefix = key.split('_')[0]
            qualifier_key = f"db_{db_prefix}"
            qualifiers[qualifier_key] = value

    return SeqFeature(location=location, type="CDS", qualifiers=qualifiers)

def generate_gbk(samples_annotations, database_list, samples_and_paths):
    os.makedirs("GBK", exist_ok=True)  # Ensure the output directory exists

    for sample, annotations in samples_annotations.items():
        print(f"Processing sample: {sample}")
        fna_file_path = samples_and_paths.get(sample)

        if not fna_file_path or not os.path.exists(fna_file_path):
            print(f"No .fna file found for sample {sample}. Skipping...")
            continue

        # Index the .fna file for efficient access
        sequence_index = SeqIO.index(fna_file_path, "fasta")

        # Initialize a SeqRecord for the sample
        sample_seq_record = SeqRecord(Seq(""), id=sample, name="", description=f"Annotations for {sample}")
        
        # Extract metadata from one of the annotations (assuming consistency)
        metadata = annotations[0] if annotations else {}
        sample_seq_record.annotations["source"] = "Your Source"
        sample_seq_record.annotations["molecule_type"] = "DNA"
        sample_seq_record.annotations["organism"] = metadata.get('taxonomy', 'Not Available')
        sample_seq_record.annotations["comment"] = f"Completeness: {metadata.get('Completeness', 'NA')}; Contamination: {metadata.get('Contamination', 'NA')}; Taxonomy: {metadata.get('taxonomy', 'Not Available')}"
        
        # Iterate through annotations, creating features for each
        for annotation in annotations:
            query_id = annotation["query_id"]
            if query_id in sequence_index:
                feature = generate_gbk_feature(annotation)
                sample_seq_record.features.append(feature)
            else:
                print(f"Warning: Query ID {query_id} not found in .fna file for sample {sample}.")

        # Determine the output file path
        output_filename = os.path.join("GBK", f"{sample}.gbk")
        
        # Write the SeqRecord to a GBK file
        with open(output_filename, "w") as output_handle:
            SeqIO.write([sample_seq_record], output_handle, "genbank")
        print(f"GBK file generated for sample {sample}: {output_filename}")



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

    # Check if GFF generation is requested and call generate_gff
    if args.gff:
        generate_gff(samples_annotations, args.database_list)

    # Check if GBK generation is requested and call generate_gbk
    if args.gbk:
        generate_gbk(samples_annotations, args.database_list, samples_and_paths)

if __name__ == "__main__":
    main()
