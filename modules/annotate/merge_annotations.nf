process MERGE_ANNOTATIONS {

    input:
    path( old_annotations, stageAs: "old_annotations.tsv" )
    path( new_annotations, stageAs: "new_annotations.tsv" )


    output:
    path "merged_combined_annotations.tsv", emit: merged_annots_out

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Function to print the first few lines of a file for debugging
    def print_file_header(file_path, num_lines=5):
        with open(file_path, 'r') as file:
            for _ in range(num_lines):
                print(file.readline())

    # Paths to your files
    old_annotations_path = "old_annotations.tsv"
    new_annotations_path = "new_annotations.tsv"

    # Print the headers of the files for debugging
    print("Old annotations file header:")
    print_file_header(old_annotations_path)

    print("New annotations file header:")
    print_file_header(new_annotations_path)

    # Function to merge or concatenate values based on conditions
    def merge_values(val1, val2):
        if pd.isnull(val1):
            return val2
        if pd.isnull(val2):
            return val1
        return val1 if val1 == val2 else f"{val1}; {val2}"

    # Load the TSV files into DataFrames
    old_df = pd.read_csv("old_annotations.tsv", sep='\t', dtype=str)
    new_df = pd.read_csv("new_annotations.tsv", sep='\t', dtype=str)

    # Perform an outer join on the two DataFrames using 'query_id' and 'sample' as keys
    merged_df = pd.merge(old_df, new_df, on=['query_id', 'sample'], how='outer', suffixes=('_old', '_new'))

    # Process columns to merge or concatenate values
    for col in old_df.columns:
        if col in new_df.columns and col not in ['query_id', 'sample']:
            merged_df[col] = merged_df.apply(lambda x: merge_values(x[col + '_old'], x[col + '_new']), axis=1)
        elif col + '_old' in merged_df.columns:
            merged_df[col] = merged_df[col + '_old']
        elif col + '_new' in merged_df.columns:
            merged_df[col] = merged_df[col + '_new']

    # Drop the _old and _new columns
    columns_to_drop = [col for col in merged_df if '_old' in col or '_new' in col]
    merged_df.drop(columns_to_drop, axis=1, inplace=True)

    # Save the merged DataFrame to a new file
    merged_file_path = "merged_combined_annotations.tsv"
    merged_df.to_csv(merged_file_path, sep='\t', index=False)


    """
}
