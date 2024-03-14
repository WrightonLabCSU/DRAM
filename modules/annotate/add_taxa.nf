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

    # Get the name of the first column in ch_taxa data
    first_column_name = ch_taxa_data.columns[0]

    # Replace "." with "-" in the first column of ch_taxa for matching
    ch_taxa_data[first_column_name] = ch_taxa_data[first_column_name].str.replace(".", "-")

    # Merge data based on the sample column
    merged_data = pd.merge(combined_annotations, ch_taxa_data[[first_column_name, 'classification']], left_on="sample", right_on=first_column_name, how="left")

    # Drop the additional first column from ch_taxa in the merged data
    merged_data.drop(columns=[first_column_name], inplace=True)

    # Rename the "classification" column to "Classification"
    merged_data.rename(columns={"classification": "taxonomy"}, inplace=True)

    # Save the updated data to annots_taxa.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
