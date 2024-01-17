def process_trnascan_output(input_file, output_file, sample_name):
    # Read the input file into a DataFrame
    trna_frame = pd.read_csv(input_file, sep="\t", skiprows=[0, 2])

    # Strip leading and trailing spaces from column names
    trna_frame.columns = trna_frame.columns.str.strip()

    # Add a new "sample" column and populate it with the sample_name value
    trna_frame.insert(0, "sample", sample_name)

    # Print column names for debugging
    print("Original column names:")
    print(trna_frame.columns)

    # Remove the "Note" column if present
    trna_frame = trna_frame.drop(columns=["Note"], errors="ignore")

    # Print column names after removing "Note" column
    print("Column names after removing 'Note':")
    print(trna_frame.columns)

    # Keep only the first occurrence of "Begin" and "End" columns
    trna_frame = trna_frame.loc[:, ~trna_frame.columns.duplicated(keep='first')]

    # Print column names after keeping the first occurrence of "Begin" and "End"
    print("Column names after keeping the first occurrence of 'Begin' and 'End':")
    print(trna_frame.columns)

    # Reorder columns
    columns_order = ["sample", "Name", "tRNA #", "begin", "end", "type", "codon", "score"]

    # Rename specified columns
    trna_frame = trna_frame.rename(columns={
        "Name": "query_id",
        "Begin": "begin",
        "End": "end",
        "Type": "type",
        "Codon": "codon",
        "Score": "score"
    })

    # Print final column names
    print("Final column names:")
    print(trna_frame.columns)

    # Write the processed DataFrame to the output file
    trna_frame.to_csv(output_file, sep="\t", index=False)