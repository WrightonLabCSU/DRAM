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
    import numpy as np

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

    # Check for existing Completeness and Contamination columns in combined_annotations
    columns_to_merge = ["Completeness", "Contamination"]
    for col in columns_to_merge:
        if col in combined_annotations.columns:
            # Check if all values are NaN in this column
            if combined_annotations[col].isnull().all():
                # Prepare to replace NaN values with checkm_data values
                continue
            else:
                # Remove the column from the list of columns to merge from checkm_data
                columns_to_merge.remove(col)

    # Adjust checkm_data to include only the necessary columns
    checkm_data_filtered = checkm_data[[first_column_name] + columns_to_merge]

    # Merge data based on the sample column
    merged_data = pd.merge(combined_annotations, checkm_data_filtered, left_on="sample", right_on=first_column_name, how="left")

    # Drop the additional first column if it's not the sample column
    if first_column_name != "sample":
        merged_data.drop(columns=[first_column_name], inplace=True)

    # Fill NaN values for merged columns if they were originally present
    if len(columns_to_merge) > 0:
        merged_data.update(checkm_data[columns_to_merge])

    # Save the updated data to annots_bin_quality.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
