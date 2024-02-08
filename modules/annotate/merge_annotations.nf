process MERGE_ANNOTATIONS {

    input:
    path( old_annotations, stageAs: "old_annotations.tsv" )
    path( new_annotations, stageAs: "new_annotations.tsv" )


    output:
    path "merged_combined_annotations.tsv", emit: merged_annots_out

    script:
    """
    #!/usr/bin/env python

    #!/usr/bin/env python

    import pandas as pd

    # Load the existing annotations file into a DataFrame
    existing_df = pd.read_csv("old_annotations.tsv", sep='\t')

    # Load the new annotations file into a DataFrame
    user_df = pd.read_csv("new_annotations.tsv", sep='\t')

    # Merge the two DataFrames based on the 'query_id' column
    merged_df = pd.merge(existing_df, user_df, on='query_id', how='outer', suffixes=('_old', '_new'))

    # Reorder columns as per the specified order
    columns_order = ['query_id', 'sample', 'start_position', 'end_position', 'strandedness']
    merged_df = merged_df.reindex(columns=columns_order + sorted(set(merged_df.columns) - set(columns_order)))

    # Save the merged DataFrame to a new file
    merged_file_path = "merged_combined_annotations.tsv"
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")

    """
}
