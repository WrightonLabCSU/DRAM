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

    def merge_annotations(old_file, new_file, output_file):
        # Load the TSV files into DataFrames
        old_df = pd.read_csv(old_file, sep='\t', dtype=str)
        new_df = pd.read_csv(new_file, sep='\t', dtype=str)

        # Merge the two DataFrames on 'query_id' and 'sample' with an outer join
        merged_df = pd.merge(old_df, new_df, on=['query_id', 'sample'], how='outer', suffixes=('_old', '_new'))

        # Iterate through the columns and merge them as necessary
        for column in merged_df.columns:
            if '_old' in column or '_new' in column:
                # Get the original column name without suffix
                original_column = column.replace('_old', '').replace('_new', '')
                # Check if the column exists in both DataFrames and merge if needed
                if (original_column + '_old') in merged_df and (original_column + '_new') in merged_df:
                    old_col = merged_df[original_column + '_old']
                    new_col = merged_df[original_column + '_new']
                    # Combine the columns if both non-null and different, otherwise take one value
                    merged_df[original_column] = old_col.where(old_col == new_col, old_col + ';' + new_col)
                    # Drop the temporary suffix columns
                    merged_df.drop([original_column + '_old', original_column + '_new'], axis=1, inplace=True)
                elif '_old' in column:
                    # Rename the column to the original name if it only exists in the old DataFrame
                    merged_df.rename(columns={column: original_column}, inplace=True)
                elif '_new' in column:
                    # Rename the column to the original name if it only exists in the new DataFrame
                    merged_df.rename(columns={column: original_column}, inplace=True)

        # Remove any duplicate rows based on 'query_id' and 'sample'
        merged_df.drop_duplicates(subset=['query_id', 'sample'], inplace=True)

        # Save the merged DataFrame to the output file
        merged_df.to_csv(output_file, sep='\t', index=False)

    # Paths to your files
    old_annotations_path = "old_annotations.tsv"
    new_annotations_path = "new_annotations.tsv"
    merged_file_path = "merged_combined_annotations.tsv"

    # Merge the annotations and save to file
    merge_annotations(old_annotations_path, new_annotations_path, merged_file_path)


    """
}
