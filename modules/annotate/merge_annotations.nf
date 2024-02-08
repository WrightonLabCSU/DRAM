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

    # Load the existing annotations file into a DataFrame
    existing_df = pd.read_csv("old_annotations.tsv", sep='\t')

    # Load the new annotations file into a DataFrame
    user_df = pd.read_csv("new_annotations.tsv", sep='\t')

    # Perform an outer join on the two DataFrames using 'query_id' and 'sample' as keys
    merged_df = pd.merge(existing_df, user_df, on=['query_id', 'sample'], how='outer', suffixes=('_old', '_new'))

    # Function to merge or concatenate values based on conditions
    def merge_values(val1, val2):
        if pd.isnull(val1):
            return val2
        if pd.isnull(val2):
            return val1
        return val1 if val1 == val2 else f"{val1}; {val2}"

    # Process columns to merge or concatenate values
    for col in existing_df.columns:
        if col in user_df.columns and col not in ['query_id', 'sample']:
            merged_df[col] = merged_df.apply(lambda x: merge_values(x[col + '_old'], x[col + '_new']), axis=1)
        elif col + '_old' in merged_df.columns:
            merged_df[col] = merged_df[col + '_old']
        elif col + '_new' in merged_df.columns:
            merged_df[col] = merged_df[col + '_new']

    # Drop the _old and _new columns
    columns_to_drop = [col for col in merged_df if '_old' in col or '_new' in col]
    merged_df.drop(columns_to_drop, axis=1, inplace=True)

    # Reorder columns as per the specified order
    columns_order = ['query_id', 'sample', 'start_position', 'end_position', 'strandedness']
    merged_df = merged_df.reindex(columns=columns_order + sorted(set(merged_df.columns) - set(columns_order)))

    # Save the merged DataFrame to a new file
    merged_file_path = "merged_combined_annotations.tsv"
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")

    """
}
