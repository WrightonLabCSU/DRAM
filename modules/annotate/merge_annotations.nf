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

    def assign_rank(row):
        rank = 'E'
        if row.get('kegg_bitScore', 0) > 350:
            rank = 'A'
        elif row.get('uniref_bitScore', 0) > 350:
            rank = 'B'
        elif row.get('kegg_bitScore', 0) > 60 or row.get('uniref_bitScore', 0) > 60:
            rank = 'C'
        elif any(row.get(f"{db}_bitScore", 0) > 60 for db in ['pfam', 'dbcan', 'merops']):
            rank = 'D'
        return rank

    def convert_bit_scores_to_numeric(df):
        for col in df.columns:
            if "_bitScore" in col:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        return df

    # Directory where the annotations TSV files are staged
    annotations_dir = "annotations"

    # Initialize an empty DataFrame for merging all annotations
    merged_df = pd.DataFrame()

    # Iterate over each TSV file in the directory
    for file_path in glob.glob(os.path.join(annotations_dir, '*.tsv')):
        # Load the current TSV file into a DataFrame
        current_df = pd.read_csv(file_path, sep='\t')
        
        # Convert bit scores to numeric for accurate comparisons and ranking
        current_df = convert_bit_scores_to_numeric(current_df)

        # If merged_df is empty, just copy the first file
        if merged_df.empty:
            merged_df = current_df
        else:
            # Perform an outer join merge with the new DataFrame
            merged_df = pd.merge(merged_df, current_df, how='outer', on=['query_id', 'sample'], suffixes=('', '_duplicate'))

    # Recalculate 'rank' for each row after merge, ensuring all columns used for ranking are numeric
    merged_df['rank'] = merged_df.apply(assign_rank, axis=1)

    # After merging, handle duplicate columns (if any) by merging their values
    for col in [col for col in merged_df.columns if '_duplicate' in col]:
        original_col = col.replace('_duplicate', '')
        # Merge values of the original and duplicate columns, then drop the duplicate
        merged_df[original_col] = merged_df.apply(lambda x: x[original_col] if pd.notnull(x[original_col]) else x[col], axis=1)
        merged_df.drop(columns=[col], inplace=True)

    # Ensure 'rank' follows 'strandedness'
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    other_columns = [col for col in merged_df.columns if col not in base_columns]
    merged_df = merged_df[base_columns + other_columns]

    # Save the merged DataFrame to a new file
    merged_file_path = "raw-merged-annotations.tsv"
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")
    """
}
