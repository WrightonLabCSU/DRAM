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

    # Load combined_annotations.tsv
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')
    # Assume genome/sample names are in the first column
    combined_annotations_first_col_name = combined_annotations.columns[0]

    # Load ch_taxa TSV
    ch_taxa_path = "${ch_taxa}"
    ch_taxa_data = pd.read_csv(ch_taxa_path, sep='\t')
    # Assume genome/sample names are in the first column of ch_taxa as well
    ch_taxa_first_col_name = ch_taxa_data.columns[0]

    # Standardize sample names in both dataframes for accurate merging
    combined_annotations[combined_annotations_first_col_name] = combined_annotations[combined_annotations_first_col_name].str.replace(".", "-", regex=False)
    ch_taxa_data[ch_taxa_first_col_name] = ch_taxa_data[ch_taxa_first_col_name].str.replace(".", "-", regex=False)

    # Merge taxonomy information based on sample names
    # Ensure only relevant columns from ch_taxa_data are merged to avoid unnecessary data duplication
    merged_data = pd.merge(combined_annotations, ch_taxa_data[[ch_taxa_first_col_name, 'classification']], left_on=combined_annotations_first_col_name, right_on=ch_taxa_first_col_name, how="left")

    # After merging, if there's an existing "taxonomy" column, we need to handle it appropriately
    if "taxonomy" in merged_data.columns:
        # If the merged 'classification' data is not null, update the 'taxonomy' column with it
        merged_data.loc[~merged_data['classification'].isnull(), 'taxonomy'] = merged_data['classification']
        # Drop the 'classification' column as its data is now merged into 'taxonomy'
        merged_data.drop(columns=['classification'], inplace=True)
    else:
        # If there's no existing "taxonomy" column, rename "classification" to "taxonomy"
        merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")
    """
}
