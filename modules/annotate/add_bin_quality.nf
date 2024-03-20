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

    # Define the column from combined_annotations that corresponds to the sample name. Adjust as necessary.
    annotations_columns = ['query_id', 'sample', 'taxonomy', 'Completeness', 'Contamination']  # Add or remove columns as needed based on your use case

    # Define columns to load from the checkm (bin quality) file
    quality_columns = ['Name', 'Completeness', 'Contamination']  # Assuming 'Name' is the column for bin/sample identifiers

    # Load combined_annotations.tsv with only necessary columns
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t', usecols=annotations_columns)

    # Load the checkm TSV (bin quality file) with only necessary columns
    checkm_path = "${ch_bin_quality}"
    checkm_data = pd.read_csv(checkm_path, sep='\t', usecols=quality_columns)

    # Standardize sample identifiers by replacing "." with "-" if necessary in both dataframes
    # Assuming the 'sample' column in combined_annotations and 'Name' in checkm_data need this adjustment
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    checkm_data["Name"] = checkm_data["Name"].str.replace(".", "-")

    # Merge operation: Proceed only if "Completeness" and "Contamination" columns are not already appropriately populated in combined_annotations
    if not combined_annotations[["Completeness", "Contamination"]].notnull().all().all():
        # Merge based on the 'sample' column in combined_annotations and 'Name' in checkm_data
        merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on="Name", how="left")
        # Drop the 'Name' column from the merge as it's redundant with 'sample'
        merged_data.drop(columns=["Name"], inplace=True)
    else:
        print("Completeness or Contamination information already exists. Skipping addition of bin quality data.")
        merged_data = combined_annotations

    # Save the updated data back to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
