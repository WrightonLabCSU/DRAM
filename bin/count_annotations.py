import pandas as pd
import click

def count_annotations(input_file, output_file):
    # Read the data from the input file
    data = pd.read_csv(input_file, sep='\t')

    # Identify columns containing database IDs and EC numbers, excluding 'query_id'
    relevant_columns = [col for col in data.columns if (col.endswith("_id") or col.endswith("_EC")) and col != 'query_id']

    # Create an empty DataFrame to store the combined data
    gene_id_data = pd.DataFrame(columns=['sample', 'gene_id'])

    # Process each relevant column, keeping entries as is without splitting
    for col in relevant_columns:
        # Keep the entries as is
        exploded_data = data[['sample', col]].copy()
        exploded_data.columns = ['sample', 'gene_id']
        
        # Append the processed data to the gene_id_data DataFrame
        gene_id_data = pd.concat([gene_id_data, exploded_data], ignore_index=True)

    # Group by "gene_id" and "sample", then count the occurrences
    table = gene_id_data.groupby(['gene_id', 'sample']).size().unstack(fill_value=0)

    # Save the resulting table to a TSV file
    table.to_csv(output_file, sep='\t')

@click.command()
@click.argument('input_file')
@click.argument('output_file')
def main(input_file, output_file):
    """Count occurrences of gene_ids and EC numbers for each sample."""
    count_annotations(input_file, output_file)


if __name__ == "__main__":
    main()
