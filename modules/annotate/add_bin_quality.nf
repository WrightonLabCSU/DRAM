process ADD_BIN_QUALITY {

    errorStrategy 'finish'

    input:
    file( combined_annotations )

    output:
    path("annots_bin_quality.tsv"), emit: annots_bin_quality_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Load combined_annotations.tsv
    combined_annotations_path = "combined_annotations.tsv"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')

    # Load checkm TSV
    checkm_path = "${params.bin_quality}"
    checkm_columns = ["Name", "Completeness", "Contamination"]
    checkm_data = pd.read_csv(checkm_path, sep='\t', usecols=checkm_columns)

    # Replace "." with "-" in the sample column for comparison
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")

    # Merge data based on the sample column
    merged_data = pd.merge(combined_annotations, checkm_data, left_on="sample", right_on="Name", how="left")

    # Drop the additional "Name" column
    merged_data.drop(columns=["Name"], inplace=True)

    # Save the updated data to annots_bin_quality.tsv
    output_path = "annots_bin_quality.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")


    """
}
