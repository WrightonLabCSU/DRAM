import pandas as pd
import argparse

def distill_summary(combined_annotations_file, genome_summary_form_file, output_file, add_modules):
    # Read input files
    combined_annotations = pd.read_csv(combined_annotations_file, sep='\t')
    genome_summary_form = pd.read_csv(genome_summary_form_file, sep='\t')

    # Process genome_summary_form
    gene_id_values = combined_annotations.filter(regex='_id$', axis=1).melt()['value'].unique()
    filtered_genome_summary = genome_summary_form[genome_summary_form['gene_id'].isin(gene_id_values)]

    # Merge combined_annotations and filtered genome_summary_form
    merged_data = pd.merge(combined_annotations, filtered_genome_summary, left_on='gene_id', right_on='gene_id', how='left')

    # Process add_modules
    for add_module_file in add_modules:
        add_module_data = pd.read_csv(add_module_file, sep='\t')
        merged_data = pd.merge(merged_data, add_module_data, left_on='gene_id', right_on='gene_id', how='left')

    # Write to output file
    merged_data.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate distill_summary.tsv')
    parser.add_argument('--combined_annotations', type=str, help='Path to combined_annotations file')
    parser.add_argument('--genome_summary_form', type=str, help='Path to genome_summary_form file')
    parser.add_argument('--output', type=str, help='Path to output file (distill_summary.tsv)')
    parser.add_argument('--add_module1', type=str, default='empty', help='Path to add_module1 file (optional)')
    parser.add_argument('--add_module2', type=str, default='empty', help='Path to add_module2 file (optional)')
    parser.add_argument('--add_module3', type=str, default='empty', help='Path to add_module3 file (optional)')
    parser.add_argument('--add_module4', type=str, default='empty', help='Path to add_module4 file (optional)')
    parser.add_argument('--add_module5', type=str, default='empty', help='Path to add_module5 file (optional)')

    args = parser.parse_args()

    add_modules = [args.add_module1, args.add_module2, args.add_module3, args.add_module4, args.add_module5]

    distill_summary(args.combined_annotations, args.genome_summary_form, args.output, add_modules)
