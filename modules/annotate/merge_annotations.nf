process MERGE_ANNOTATIONS {

    input:
    path old_annotations
    path new_annotations

    output:
    path "merged_combined_annotations.tsv", emit: merged_annots_out

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    # Load the existing annotations file into a DataFrame
    existing_df = pd.read_csv(${old_annotations}, sep='\\t')

    # Load the new annotations file into a DataFrame
    user_df = pd.read_csv(${new_annotations}, sep='\\t')

    # Merge the two DataFrames based on the 'query_id' column
    merged_df = pd.merge(existing_df, user_df, on='query_id', how='outer')

    # Save the merged DataFrame to a new file
    merged_file_path = "${merged_annots_out}"
    merged_df.to_csv(merged_file_path, sep='\\t', index=False)

    print(f"Merged annotations saved to {merged_file_path}")
    """
}
