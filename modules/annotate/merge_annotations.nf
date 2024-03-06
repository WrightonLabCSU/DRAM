process MERGE_ANNOTATIONS {

    input:
    path( ch_annotations, stageAs: "annotations/*" )



    output:
    path "raw-merged-annotations.tsv", emit: ch_merged_annots_out

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import os
    import glob

    # Directory where the annotations TSV files are staged
    annotations_dir = "annotations"

    # Initialize an empty DataFrame for merging all annotations
    merged_df = pd.DataFrame()

    # Iterate over each TSV file in the directory
    for file_path in glob.glob(os.path.join(annotations_dir, '*.tsv')):
        # Load the current TSV file into a DataFrame
        current_df = pd.read_csv(file_path, sep='\t')
        
        # If merged_df is empty, just copy the first file
        if merged_df.empty:
            merged_df = current_df
        else:
            # Perform an outer join merge with the new DataFrame
            merged_df = pd.merge(merged_df, current_df, how='outer', on=['query_id', 'sample'], suffixes=('', '_duplicate'))

    # After merging, handle duplicate columns (if any) by merging their values
    for col in [col for col in merged_df.columns if '_duplicate' in col]:
        original_col = col.replace('_duplicate', '')
        # Merge values of the original and duplicate columns, then drop the duplicate
        merged_df[original_col] = merged_df.apply(lambda x: x[original_col] if pd.notnull(x[original_col]) else x[col], axis=1)
        merged_df.drop(columns=[col], inplace=True)

    # Save the merged DataFrame to a new file
    merged_file_path = "raw-merged-annotations.tsv"
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")

    """
}
