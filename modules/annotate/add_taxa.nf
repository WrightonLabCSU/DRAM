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

    # Define necessary columns for loading. Adjust column names as needed.
    annotations_columns = ['sample', 'other_columns_if_necessary']
    taxa_columns = ['sample_name', 'classification'] # Adjust 'sample_name' to the actual name of the first column in ch_taxa

    # Load combined_annotations.tsv with only necessary columns
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t', usecols=annotations_columns)

    # Load ch_taxa TSV with only necessary columns
    ch_taxa_path = "${ch_taxa}"
    ch_taxa_data = pd.read_csv(ch_taxa_path, sep='\t', usecols=taxa_columns)

    # Efficiently replace "." with "-" in the 'sample' columns
    # Assuming 'sample' is the name in both dataframes, adjust if it's different in ch_taxa
    combined_annotations['sample'] = combined_annotations['sample'].str.replace('.', '-')
    ch_taxa_data['sample_name'] = ch_taxa_data['sample_name'].str.replace('.', '-')

    # Check if 'taxonomy' column exists and has at least one non-empty value
    if 'taxonomy' in combined_annotations.columns and not combined_annotations['taxonomy'].isnull().all():
        print("Taxonomy information already exists. Skipping addition of taxonomy data.")
    else:
        # Proceed with adding taxonomy information

        # Merge data based on the sample column. Adjust 'sample_name' to match the actual column name in ch_taxa
        merged_data = pd.merge(combined_annotations, ch_taxa_data, left_on="sample", right_on="sample_name", how="left")

        # Drop the additional 'sample_name' column from ch_taxa in the merged data
        merged_data.drop(columns=['sample_name'], inplace=True)

        # Rename the "classification" column to "taxonomy"
        merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)

        # Save the updated data to raw-annotations.tsv
        output_path = "raw-annotations.tsv"
        merged_data.to_csv(output_path, sep='\t', index=False)

        print(f"Updated annotations saved to {output_path}")
    """
}
