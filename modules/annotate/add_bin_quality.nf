process ADD_BIN_QUALITY {

    errorStrategy 'finish'

    input:
    file( combined_annotations )
    file( ch_bin_quality )

    output:
    path("raw-annotations.tsv"), emit: annots_bin_quality_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Load combined_annotations.tsv
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')

    # Load checkm TSV
    checkm_path = "${ch_bin_quality}"
    checkm_data = pd.read_csv(checkm_path, sep='\t')

    # Identify the first column name dynamically
    first_column_name = checkm_data.columns[0]

    # Extract relevant columns from checkm_data
    checkm_columns = [first_column_name, "Completeness", "Contamination"]
    checkm_data = checkm_data[checkm_columns]

    # Replace "." with "-" in the sample column for comparison
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")

    # Replace "." with "-" in the first column of checkm TSV
    checkm_data[first_column_name] = checkm_data[first_column_name].str.replace(".", "-")

    # Check if "Completeness" and "Contamination" columns exist and have non-null values in combined_annotations
    if "Completeness" in combined_annotations.columns and "Contamination" in combined_annotations.columns:
        # Check if all values are NaN in these columns
        if combined_annotations["Completeness"].isnull().all() and combined_annotations["Contamination"].isnull().all():
            # If both columns are empty, prepare to replace them with checkm_data values
            completeness_contamination_to_merge = True
        else:
            # If any of the columns have non-null values, do not merge from checkm_data
            completeness_contamination_to_merge = False
    else:
        # If either or both columns don't exist, prepare to add them
        completeness_contamination_to_merge = True

    if completeness_contamination_to_merge:
        # Merge data based on the sample column
        merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on=first_column_name, how="left")
        # Drop the additional first column if it's not the sample column
        if first_column_name != "sample":
            merged_data.drop(columns=[first_column_name], inplace=True)
    else:
        # If Completeness and Contamination should not be merged, use combined_annotations directly
        merged_data = combined_annotations

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
