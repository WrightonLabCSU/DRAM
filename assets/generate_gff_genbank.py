import argparse
import csv
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Generate GFF and/or GBK files from raw annotations.")
    parser.add_argument("--gff", action='store_true', help="Generate GFF file")
    parser.add_argument("--gbk", action='store_true', help="Generate GBK file")
    parser.add_argument("--annotations", required=True, help="Path to the raw annotations file")
    return parser.parse_args()

def sanitize_description(description):
    """Replace semicolons in descriptions to avoid parsing issues."""
    return description.replace(';', ',')

def format_attributes(annotation, prioritized_keys=['kegg_id', 'kofam_id']):
    """Format and order database-specific annotations for the GFF attributes column."""
    attributes = []
    # Ensure prioritized keys are handled first if they exist
    for key in prioritized_keys:
        if key in annotation and annotation[key]:
            desc_key = f"{key}_description"
            desc = sanitize_description(annotation.get(desc_key, "NA"))
            attributes.append(f"{key}={annotation[key]};description={desc}")

    # Handle other keys in alphabetical order
    sorted_keys = sorted(set(annotation.keys()) - set(prioritized_keys) - {'query_id', 'sample', 'start_position', 'end_position', 'strandedness', 'Completeness', 'Contamination', 'taxonomy'})
    for key in sorted_keys:
        if key.endswith('_id'):
            desc_key = key.replace('_id', '_description')
            ec_key = key.replace('_id', '_EC')
            desc = sanitize_description(annotation.get(desc_key, "NA"))
            if ec_key in annotation and annotation[ec_key]:
                desc += f"; EC:{annotation[ec_key]}"
            attributes.append(f"{key}={annotation[key]};description={desc}")

    return "; ".join(attributes)

def generate_gff(samples_annotations):
    """Generate GFF files for each sample."""
    for sample, annotations in samples_annotations.items():
        with open(f"{sample}.gff", "w") as gff_file:
            # Write metadata as comments
            metadata = annotations[0]  # Assuming all entries for a sample share the same metadata
            gff_file.write(f"##gff-version 3\n")
            gff_file.write(f"# Completeness: {metadata['Completeness']}\n")
            gff_file.write(f"# Contamination: {metadata['Contamination']}\n")
            gff_file.write(f"# Taxonomy: {metadata['taxonomy']}\n")
            
            for annotation in annotations:
                # Prepare attributes
                attributes_str = format_attributes(annotation)
                # Format the GFF line
                gff_line = f"{annotation['query_id']}\t.\tgene\t{annotation['start_position']}\t{annotation['end_position']}\t.\t{'+' if annotation['strandedness'] == '+1' else '-'}\t.\t{attributes_str}\n"
                gff_file.write(gff_line)

def generate_gbk(samples_annotations):
    """Generate GBK files for each sample."""
    # Placeholder for the GBK generation code
    print("GBK generation not implemented yet.")

def main():
    args = parse_arguments()

    # Load annotations and organize by sample
    samples_annotations = defaultdict(list)
    with open(args.annotations, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            samples_annotations[row['sample']].append(row)

    # Generate GFF and/or GBK files based on flags
    if args.gff:
        generate_gff(samples_annotations)
    if args.gbk:
        generate_gbk(samples_annotations)

if __name__ == "__main__":
    main()
