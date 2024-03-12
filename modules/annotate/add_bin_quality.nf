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

    # Ensure "Completeness" and "Contamination" information is updated only if missing in combined_annotations
    # Merge data based on the sample column, prioritizing existing data in combined_annotations
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on=first_column_name, how="left")

    # For each column in checkm_columns (excluding the matching column), update only if NaN in combined_annotations
    for col in ["Completeness", "Contamination"]:
        # Use combined_annotations' values if present, otherwise fill in from checkm_data
        if col in combined_annotations.columns:
            merged_data[col] = merged_data[col + "_x"].combine_first(merged_data[col + "_y"])
        else:  # If the column doesn't exist in combined_annotations, simply use checkm_data's values
            merged_data[col] = merged_data[col + "_y"]

        # Cleanup temporary merge columns
        merged_data.drop(columns=[col + "_x", col + "_y"], inplace=True, errors='ignore')

    # Ensure the first column from checkm_data is dropped if it was merged and is no longer needed
    if first_column_name != "sample":
        merged_data.drop(columns=[first_column_name], inplace=True, errors='ignore')

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")


    """
}
