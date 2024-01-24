process ADD_TAXA {

    errorStrategy 'finish'

    input:
    file( combined_annotations )
    file( ch_taxa )

    output:
    path("annots_taxa.tsv"), emit: annots_taxa_out, optional: true

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

    # Rename the "classification" column to "Classification"
    merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)

    # Save the updated data to annots_taxa.tsv
    output_path = "annots_taxa.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
