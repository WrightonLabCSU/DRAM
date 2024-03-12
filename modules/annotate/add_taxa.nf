process ADD_TAXA {

    errorStrategy 'finish'

    input:
    file( combined_annotations )
    file( ch_taxa )

    output:
    path("raw-annotations.tsv"), emit: annots_taxa_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Load combined_annotations.tsv
    combined_annotations_path = "${combined_annotations}"
    combined_annotations = pd.read_csv(combined_annotations_path, sep='\t')

    # Load ch_taxa TSV
    ch_taxa_path = "${ch_taxa}"
    ch_taxa_data = pd.read_csv(ch_taxa_path, sep='\t')

    # Replace "." with "-" in the sample column for comparison
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")

    # Replace "." with "-" in the "user_genome" column of ch_taxa for matching
    ch_taxa_data["user_genome"] = ch_taxa_data["user_genome"].str.replace(".", "-")

    # Check if "taxonomy" column exists and has non-null values
    if "taxonomy" in combined_annotations.columns:
        if not combined_annotations["taxonomy"].isnull().all():
            # If "taxonomy" exists and has non-null values, don't merge taxonomy from ch_taxa_data
            taxonomy_to_merge = False
        else:
            # If "taxonomy" exists but only contains null values, prepare to replace those with ch_taxa_data
            taxonomy_to_merge = True
    else:
        # If "taxonomy" doesn't exist, add it
        taxonomy_to_merge = True

    if taxonomy_to_merge:
        # Merge data based on the sample column
        merged_data = pd.merge(combined_annotations, ch_taxa_data[['user_genome', 'classification']], left_on="sample", right_on="user_genome", how="left")
        # Drop the additional "user_genome" column
        merged_data.drop(columns=["user_genome"], inplace=True)
        # Rename the "classification" column to "taxonomy"
        merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)
    else:
        # If taxonomy should not be merged, use combined_annotations directly
        merged_data = combined_annotations

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
