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

    # Load the entire combined_annotations.tsv
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')

    # Define columns to load from the ch_bin_quality (checkm) file
    # Assuming 'Name' is the column for bin/sample identifiers in ch_bin_quality
    quality_columns = ['Name', 'Completeness', 'Contamination']

    # Load the ch_bin_quality TSV with only necessary columns
    checkm_path = "${ch_bin_quality}"
    checkm_data = pd.read_csv(checkm_path, sep='\t', usecols=quality_columns)

    # Standardize sample identifiers by replacing "." with "-" if necessary in both dataframes
    # Assuming the 'sample' column in combined_annotations corresponds to 'Name' in checkm_data
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    checkm_data["Name"] = checkm_data["Name"].str.replace(".", "-")

    # Merge operation: Add "Completeness" and "Contamination" from checkm_data to combined_annotations
    # Merge based on the 'sample' column in combined_annotations and 'Name' in checkm_data
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on="Name", how="left")

    # Drop the 'Name' column from the merge as it's redundant with 'sample'
    merged_data.drop(columns=["Name"], inplace=True)

    # Save the updated data back to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")


    """
}
