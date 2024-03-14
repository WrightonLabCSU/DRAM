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

    # Replace "." with "-" in the first column for comparison in combined_annotations DataFrame
    combined_annotations[combined_annotations_first_col_name] = combined_annotations[combined_annotations_first_col_name].str.replace(".", "-", regex=False)

    # Since taxonomy information is always in a column named 'classification', no need to check column name in ch_taxa_data
    # However, we need to ensure that we replace "." with "-" in the first column for comparison, if needed
    ch_taxa_first_col_name = ch_taxa_data.columns[0]
    ch_taxa_data[ch_taxa_first_col_name] = ch_taxa_data[ch_taxa_first_col_name].str.replace(".", "-", regex=False)

    # Check if "taxonomy" column exists and has non-null values in combined_annotations
    if "taxonomy" in combined_annotations.columns and not combined_annotations["taxonomy"].isnull().all():
        taxonomy_to_merge = False
    else:
        taxonomy_to_merge = True

    if taxonomy_to_merge:
        # Merge data based on the first column which represents the genome/sample name
        merged_data = pd.merge(combined_annotations, ch_taxa_data[[ch_taxa_first_col_name, 'classification']], left_on=combined_annotations_first_col_name, right_on=ch_taxa_first_col_name, how="left")
        # Drop the additional first column from ch_taxa_data that was used for the merge
        merged_data.drop(columns=[ch_taxa_first_col_name], inplace=True)
        # Rename the "classification" column to "taxonomy"
        merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)
    else:
        # If taxonomy should not be merged, use combined_annotations directly
        merged_data = combined_annotations

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")
    """
}
