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
    existing_df = pd.read_csv(existing_file_path, sep='	')

    # Load the user-provided combined_annotations.tsv file into a DataFrame
    user_df = pd.read_csv(user_file_path, sep='	')

    # Merge the two DataFrames based on the 'query_id' column
    merged_df = pd.merge(existing_df, user_df, on='query_id', how='outer', suffixes=('_existing', '_user'))

    # Initialize a list to keep track of columns processed
    processed_cols = []

    # Iterate over columns that exist in both files
    for col in merged_df.columns:
        # Check if the column is present in both files and if it's not 'query_id'
        if col in existing_df.columns and col in user_df.columns and col != 'query_id':
            merged_df[col] = merged_df.apply(
                lambda row: f"{row[col + '_existing']}; {row[col + '_user']}" if row[col + '_existing'] != row[col + '_user'] else row[col + '_existing'],
                axis=1
            )
            # Drop the intermediate columns '_existing' and '_user'
            merged_df = merged_df.drop([col + '_existing', col + '_user'], axis=1)
            processed_cols.extend([col + '_existing', col + '_user'])

    # Remove the 'query_id' columns from the processed list
    processed_cols = [col for col in processed_cols if col != 'query_id_existing' and col != 'query_id_user']

    # Remove duplicate 'query_id' column if it exists
    if 'query_id_user' in merged_df.columns and 'query_id_existing' in merged_df.columns:
        merged_df = merged_df.drop('query_id_user', axis=1)

    # Rename the 'query_id_existing' column to 'query_id'
    if 'query_id_existing' in merged_df.columns:
        merged_df = merged_df.rename(columns={'query_id_existing': 'query_id'})

    # Save the merged DataFrame to a new file
    merged_file_path = 'merged_combined_annotations.tsv'
    merged_df.to_csv(merged_file_path, sep='	', index=False)

    print(f"Merged annotations saved to {merged_file_path}")

    """
}
