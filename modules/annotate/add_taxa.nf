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
    ch_taxa_data = pd.read_csv(ch_taxa_path, sep='\t', dtype=str)  # Ensure all data is read as string to avoid dtype issues

    # Replace "." with "-" for consistency
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    ch_taxa_data["user_genome"] = ch_taxa_data["user_genome"].str.replace(".", "-")

    # Merge taxonomy, ensuring 'rank' is preserved
    merged_data = pd.merge(combined_annotations, ch_taxa_data[['user_genome', 'classification']], left_on="sample", right_on="user_genome", how="left", suffixes=('', '_taxa'))

    # Update 'taxonomy' with 'classification' if it's NaN in the original data
    if 'taxonomy' in merged_data and 'classification' in merged_data:
        merged_data['taxonomy'] = merged_data['taxonomy'].fillna(merged_data['classification'])
        merged_data.drop(['classification', 'user_genome'], axis=1, inplace=True, errors='ignore')

    # Save the updated data
    merged_data.to_csv("raw-annotations.tsv", sep='\t', index=False)

    """
}
