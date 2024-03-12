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
    checkm_data = pd.read_csv(checkm_path, sep='\t', dtype=str)  # Ensure all data is read as string to avoid dtype issues

    # Replace "." with "-" in the sample column for consistent naming
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    checkm_data[checkm_data.columns[0]] = checkm_data[checkm_data.columns[0]].str.replace(".", "-")

    # Merge while keeping the 'rank' column intact
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on=checkm_data.columns[0], how="left")

    # Fill in "Completeness" and "Contamination" from checkm_data if they are NaN in combined_annotations
    for col in ["Completeness", "Contamination"]:
        if col in merged_data.columns:
            merged_data[col] = merged_data[col].fillna(merged_data[f'{col}_y'])
            merged_data.drop([f'{col}_y'], axis=1, inplace=True, errors='ignore')

    # Clean up any temporary columns and ensure 'rank' column is not affected
    merged_data.drop(list(merged_data.filter(regex='_x$')), axis=1, inplace=True, errors='ignore')
    merged_data.drop(list(merged_data.filter(regex='_y$')), axis=1, inplace=True, errors='ignore')

    # Save the updated data
    merged_data.to_csv("raw-annotations.tsv", sep='\t', index=False)

    """
}
