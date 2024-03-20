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

    # Define necessary columns for loading
    annotations_columns = ['sample', 'other_relevant_columns']  # Update other_relevant_columns as needed
    quality_columns = ['Name', 'Completeness', 'Contamination']  # Assuming 'Name' is the sample identifier in ch_bin_quality

    # Load combined_annotations.tsv with only necessary columns
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\\t', usecols=annotations_columns)

    # Load checkm TSV with only necessary columns
    checkm_path = "${ch_bin_quality}"
    checkm_data = pd.read_csv(checkm_path, sep='\\t', usecols=quality_columns)

    # Standardize sample identifiers by replacing "." with "-" if necessary
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    checkm_data["Name"] = checkm_data["Name"].str.replace(".", "-")

    # Determine if "Completeness" and "Contamination" need to be merged
    if "Completeness" in combined_annotations.columns and "Contamination" in combined_annotations.columns:
        if not combined_annotations["Completeness"].isnull().all() or not combined_annotations["Contamination"].isnull().all():
            print("Completeness or Contamination information already exists. Skipping addition of bin quality data.")
            merged_data = combined_annotations
        else:
            # Merge and update Completeness and Contamination
            merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on="Name", how="left")
            merged_data.drop(columns=["Name"], inplace=True)
    else:
        # Columns don't exist, proceed to merge
        merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on="Name", how="left")
        merged_data.drop(columns=["Name"], inplace=True)

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\\t', index=False)

    print(f"Updated annotations saved to {output_path}")
    """
}
