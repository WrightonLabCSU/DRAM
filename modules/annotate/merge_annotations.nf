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

    # Initialize a list to collect DataFrames
    dfs = []

    # Iterate over each TSV file in the directory
    for file_path in glob.glob(os.path.join(annotations_dir, '*.tsv')):
        # Load the current TSV file into a DataFrame
        current_df = pd.read_csv(file_path, sep='\t')
        
        # Convert bit scores to numeric for accurate comparisons and ranking
        current_df = convert_bit_scores_to_numeric(current_df)

        # Add DataFrame to the list
        dfs.append(current_df)

    # Concatenate all DataFrames in the list
    merged_df = pd.concat(dfs, ignore_index=True)

    # Recalculate 'rank' for each row after concatenation, ensuring all columns used for ranking are numeric
    merged_df['rank'] = merged_df.apply(assign_rank, axis=1)

    # Insert 'rank' after 'strandedness' and organize columns as required
    base_columns = ['query_id', 'sample', 'start_position', 'stop_position', 'strandedness', 'rank']
    other_columns = [col for col in merged_df.columns if col not in base_columns]
    merged_df = merged_df[base_columns + other_columns]

    # Save the merged DataFrame to a new file
    merged_file_path = "raw-merged-annotations.tsv"
    merged_df.to_csv(merged_file_path, sep='\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")
    """
}
