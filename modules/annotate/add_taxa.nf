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

    # Prepare for comparison adjustments
    combined_annotations["sample"] = combined_annotations["sample"].str.replace(".", "-")
    ch_taxa_data["user_genome"] = ch_taxa_data["user_genome"].str.replace(".", "-")

    # Merge taxonomy data if "taxonomy" column doesn't exist or is entirely null
    if "taxonomy" not in combined_annotations.columns or combined_annotations["taxonomy"].isnull().all():
        merged_data = pd.merge(combined_annotations, ch_taxa_data[['user_genome', 'classification']], left_on="sample", right_on="user_genome", how="left", suffixes=('', '_taxa'))
        merged_data.drop(columns=['user_genome'], inplace=True, errors='ignore')
        merged_data.rename(columns={'classification': 'taxonomy'}, inplace=True)
    else:
        merged_data = combined_annotations.copy()

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")

    """
}
