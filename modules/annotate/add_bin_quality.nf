process ADD_BIN_QUALITY {

    errorStrategy 'finish'

    input:
    file( combined_annotations )
    file( ch_bin_quality )

    output:
    path("annots_bin_quality.tsv"), emit: annots_bin_quality_out, optional: true

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

    # Merge data based on the sample column
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on=first_column_name, how="left")

    # Drop the additional first column
    merged_data.drop(columns=[first_column_name], inplace=True)

    # Save the updated data to annots_bin_quality.tsv
    output_path = "annots_bin_quality.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
