process MERGE_ANNOTATIONS {

    input:
    path( old_annotations, stageAs: "old_annotations.tsv" )
    path( new_annotations, stageAs: "new_annotations.tsv" )


    output:
    path "raw-annotations.tsv", emit: merged_annots_out

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Define the paths to the existing combined_annotations.tsv file and the user-provided file
    existing_file_path = 'old_annotations.tsv'
    user_file_path = 'new_annotations.tsv'

    # Load the existing combined_annotations.tsv file into a DataFrame
    existing_df = pd.read_csv(existing_file_path, sep='\t')

    # Load the user-provided combined_annotations.tsv file into a DataFrame
    user_df = pd.read_csv(user_file_path, sep='\t')

    # Merge the two DataFrames based on the 'query_id' column, retaining all columns from both files
    merged_df = pd.merge(existing_df, user_df, on='query_id', how='outer')

    # Iterate over columns that need special handling (i.e., present in both files)
    for col in merged_df.columns:
        # Check if the column is present in both files and if the data differs for a given 'query_id'
        if col in existing_df.columns and col in user_df.columns:
            merged_df[col] = merged_df.apply(
                lambda row: f"{row[col + '_x']}; {row[col + '_y']}" if row[col + '_x'] != row[col + '_y'] else row[col + '_x'],
                axis=1
            )
            # Drop the intermediate columns '_x' and '_y'
            merged_df = merged_df.drop([col + '_x', col + '_y'], axis=1)

    # Save the merged DataFrame to a new file
    merged_file_path = 'merged_combined_annotations.tsv'
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")


    """
}
