process MERGE_ANNOTATIONS {

    input:
    file( old_annotations, stageAs: "old_annotations.tsv" )
    file( new_annotations, stageAs: "new_annotations.tsv" )


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

    # Merge the two DataFrames based on the 'query_id' column, retaining all rows from both files
    merged_df = pd.concat([existing_df, user_df], ignore_index=True, sort=False)

    # Iterate over columns that need special handling (e.g., 'dbcan_id')
    columns_to_concatenate = ['dbcan_id']  # Add more columns as needed
    for col in columns_to_concatenate:
        merged_df[col] = merged_df.groupby('query_id')[col].transform(lambda x: '; '.join(x.unique()))

    # Save the merged DataFrame to a new file
    merged_file_path = 'merged_combined_annotations.tsv'
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")

    """
}
