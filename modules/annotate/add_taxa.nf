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
    # Assume genome/sample names are in the first column
    combined_annotations_first_col_name = combined_annotations.columns[0]

    # Load ch_taxa TSV
    ch_taxa_path = "${ch_taxa}"
    ch_taxa_data = pd.read_csv(ch_taxa_path, sep='\t')
    # Standardize sample names in ch_taxa_data for accurate merging
    ch_taxa_data.iloc[:, 0] = ch_taxa_data.iloc[:, 0].str.replace(".", "-", regex=False)

    # Merge ch_taxa_data into combined_annotations on sample names
    # Ensure we're merging on the first column of both DataFrames
    merged_data = pd.merge(combined_annotations, ch_taxa_data, left_on=combined_annotations_first_col_name, right_on=ch_taxa_data.columns[0], how="left")

    # Rename 'classification' column to 'taxonomy' if it doesn't exist or update existing 'taxonomy' column
    if 'classification' in merged_data.columns:
        if 'taxonomy' not in combined_annotations.columns:
            merged_data.rename(columns={'classification': 'taxonomy'}, inplace=True)
        else:
            # Update 'taxonomy' column with 'classification' values where 'classification' is not null
            merged_data['taxonomy'].update(merged_data.pop('classification'))

    # Remove any duplicate columns that might have resulted from the merge
    merged_data = merged_data.loc[:,~merged_data.columns.duplicated()]

    # Save the updated data to raw-annotations.tsv
    output_path = "raw-annotations.tsv"
    merged_data.to_csv(output_path, sep='\t', index=False)

    print(f"Updated annotations saved to {output_path}")
    """
}
