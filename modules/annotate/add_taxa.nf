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

    # Merge data based on the sample column
    merged_data = pd.merge(combined_annotations, ch_taxa_data[['user_genome', 'classification']], left_on="sample", right_on="user_genome", how="left")
    # Drop the additional "user_genome" column
    merged_data.drop(columns=["user_genome"], inplace=True)
    # Rename the "classification" column to "taxonomy"
    merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)

    # Check if "taxonomy" column exists
    if "taxonomy" in combined_annotations.columns:
        # Replace NaN values in the existing taxonomy column with corresponding values from ch_taxa_data
        merged_data['taxonomy'] = merged_data['taxonomy'].combine_first(ch_taxa_data.set_index('user_genome')['classification'])

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
