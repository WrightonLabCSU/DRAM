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

    # Load the ch_bin_quality TSV header to check for 'Name' or 'Bin Id'
    checkm_path = "${ch_bin_quality}"
    checkm_header = pd.read_csv(checkm_path, sep='\t', nrows=0)  # Load only the header

    # Determine the correct column name ('Name' or 'Bin Id') for bin/sample identifiers
    if 'Name' in checkm_header.columns:
        id_column = 'Name'
    elif 'Bin Id' in checkm_header.columns:
        id_column = 'Bin Id'
    else:
        raise ValueError("Input ch_bin_quality file must contain either 'Name' or 'Bin Id' column.")

    # Define columns to load from the ch_bin_quality file, including the determined id_column
    quality_columns = [id_column, 'Completeness', 'Contamination']

    # Load the ch_bin_quality TSV with only necessary columns
    checkm_data = pd.read_csv(checkm_path, sep='\t', usecols=quality_columns)

    # Standardize sample identifiers by replacing "." with "-" if necessary in both dataframes
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    checkm_data[id_column] = checkm_data[id_column].str.replace(".", "-")

    # Merge operation: Add "Completeness" and "Contamination" from checkm_data to combined_annotations
    # Merge based on the 'sample' column in combined_annotations and id_column in checkm_data
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on=id_column, how="left")

    # Drop the id_column from the merge as it's redundant with 'sample'
    merged_data.drop(columns=[id_column], inplace=True)

    # Save the updated data back to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")
    """
}
