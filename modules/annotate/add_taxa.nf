process ADD_TAXA {

    errorStrategy 'finish'

    input:
    file( combined_annotations )
    file( ch_taxa )

    output:
    path("raw-annotations.tsv"), emit: annots_taxa_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Adjust these column names based on the actual columns you need to use
    annotations_columns = ['query_id', 'sample']  # 'query_id' is included as an example; adjust as needed
    taxa_columns = ['user_genome', 'classification']

    # Corrected file paths
    combined_annotations_path = "raw-annotations.tsv"
    ch_taxa_path = "gtdbtk.bac120.summary.tsv"

    # Load combined_annotations.tsv with only necessary columns
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t', usecols=annotations_columns)

    # Load ch_taxa TSV with only necessary columns
    ch_taxa_data = pd.read_csv(ch_taxa_path, sep='\t', usecols=taxa_columns)

    # Correctly replace "." with "-" in the 'sample' column and 'user_genome' in ch_taxa
    combined_annotations['sample'] = combined_annotations['sample'].str.replace('.', '-')
    ch_taxa_data['user_genome'] = ch_taxa_data['user_genome'].str.replace('.', '-')

    # Check if 'taxonomy' column exists and has at least one non-empty value
    if 'taxonomy' in combined_annotations.columns and not combined_annotations['taxonomy'].isnull().all():
        print("Taxonomy information already exists. Skipping addition of taxonomy data.")
    else:
        # Proceed with adding taxonomy information by merging based on the corrected columns
        merged_data = pd.merge(combined_annotations, ch_taxa_data, left_on="sample", right_on="user_genome", how="left")

        # Drop the additional 'user_genome' column from ch_taxa in the merged data
        merged_data.drop(columns=['user_genome'], inplace=True)

        # Rename the "classification" column to "taxonomy"
        merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)

        # Save the updated data to raw-annotations.tsv
        output_path = "raw-annotations.tsv"  # Changed to avoid overwriting the original file
        merged_data.to_csv(output_path, sep='\t', index=False)

        print(f"Updated annotations saved to {output_path}")

    """
}
