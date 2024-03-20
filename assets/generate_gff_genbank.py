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
            
            # Metadata: Completeness, Contamination, and sanitized Taxonomy
            completeness = annotations[0].get('Completeness', 'NA')
            contamination = annotations[0].get('Contamination', 'NA')
            taxonomy = escape_gff3_value(annotations[0].get('taxonomy', 'NA').replace(';', ','))  # Sanitize taxonomy
            
            gff_file.write(f"# Completeness: {completeness}\n")
            gff_file.write(f"# Contamination: {contamination}\n")
            gff_file.write(f"# Taxonomy: {taxonomy}\n")
            
            for annotation in annotations:
                seqid = escape_gff3_value(annotation['query_id'])
                source = '.'  # Use a real source if available
                type = "gene"  # Adjust based on actual data
                start = annotation['start_position']
                end = annotation['stop_position']
                score = '.'  # Use actual score if available
                strand = '+' if annotation['strandedness'] == '+1' else '-'
                phase = '.'  # Adjust for features like CDS
                attributes_str = format_attributes(annotation, database_list)
                
                gff_line = f"{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes_str}\n"
                gff_file.write(gff_line)

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

def format_qualifiers(annotation):
    """
    Formats qualifiers for a GenBank feature, ensuring correct syntax and inclusion of essential information.
    """
    qualifiers = []
    for key, value in annotation.items():
        if key in ['gene', 'product', 'note', 'function']:  # Example qualifier keys
            value = value.replace('"', '""')  # Escape double quotes in value
            qualifiers.append(f"/{key}=\"{value}\"")
    # Add additional qualifiers as needed, using controlled vocabularies or specific formats
    return qualifiers

def format_qualifiers_gbk(qualifiers_dict):
    """
    Format qualifiers for a GenBank feature.
    """
    qualifiers_list = []
    for key, value in qualifiers_dict.items():
        # Standardize the representation of certain key qualifiers for GenBank format
        if key in ['gene', 'product', 'note', 'function']:
            qualifiers_list.append(f"/{key}=\"{value}\"")
    return qualifiers_list

def generate_gbk_feature(feature, seq_record):
    """
    Generate a SeqFeature for a GenBank file based on annotation data.
    Assumes 1-based indexing in the feature dictionary and converts to 0-based indexing for Biopython.
    """
    start = int(feature['start_position']) - 1  # Convert start to 0-based indexing
    end = int(feature['stop_position'])  # End is inclusive, no conversion needed
    strand = 1 if feature['strandedness'] == '+1' else -1  # Convert strandedness to Biopython strand
    feature_location = FeatureLocation(start, end, strand=strand)
    qualifiers = {
        'locus_tag': feature.get('query_id'),
        'product': feature.get('product', 'unknown product')
    }
    return SeqFeature(location=feature_location, type="CDS", qualifiers=qualifiers)

def generate_gbk(samples_annotations, database_list, samples_and_paths):
    print("Starting GBK generation...")
    os.makedirs("GBK", exist_ok=True)

    for sample, annotations in samples_annotations.items():
        # Retrieve the path to the sample's .fna file
        fna_file_path = samples_and_paths.get(sample)
        print(f"\nProcessing sample: {sample}")
        
        if fna_file_path and os.path.exists(fna_file_path):
            sample_sequences = SeqIO.index(fna_file_path, "fasta")  # Index the sequences for efficient access

            for annotation in annotations:
                query_id = annotation['query_id']
                if query_id in sample_sequences:
                    seq_record = SeqRecord(sample_sequences[query_id].seq, id=query_id, description="")
                    seq_record.annotations["source"] = "DRAM2"
                    seq_record.annotations["molecule_type"] = "DNA"

                    # Optional: Extract taxonomy, completeness, and contamination from annotations if applicable
                    # seq_record.annotations["organism"] = annotation.get('taxonomy', 'Not Available')
                    # Custom metadata example
                    # custom_metadata = f"Completeness: {annotation.get('Completeness', 'NA')}; Contamination: {annotation.get('Contamination', 'NA')}"
                    # seq_record.annotations["note"] = custom_metadata

                    feature_location = FeatureLocation(start=int(annotation['start_position']) - 1, end=int(annotation['stop_position']), strand=1 if annotation['strandedness'] == '+1' else -1)
                    qualifiers = format_qualifiers(annotation, database_list)
                    feature = SeqFeature(feature_location, type="gene", qualifiers=qualifiers)
                    seq_record.features.append(feature)

                    output_filename = os.path.join("GBK", f"{query_id}.gbk")
                    with open(output_filename, "w") as output_handle:
                        SeqIO.write([seq_record], output_handle, "genbank")
                    print(f"GBK file generated for query ID {query_id}: {output_filename}")
                else:
                    print(f"Sequence for query ID {query_id} not found in provided .fna file for sample {sample}.")
        else:
            print(f".fna file path not found or does not exist for sample {sample}. Skipping...")


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
