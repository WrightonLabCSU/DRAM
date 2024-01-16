process TRNA_SCAN {

    errorStrategy 'finish'

    tag { sample }

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_processed_trnas.tsv"), emit: trna_scan_out, optional: true

    script:

    script:

    """
    python3 - <<EOF
    import pandas as pd
    import subprocess
    from nextflow import Nextflow

    # Access the threads parameter using the nextflow object
    threads = Nextflow().params.threads

    # Run tRNAscan-SE
    subprocess.run(['tRNAscan-SE', '-G', '--thread', str(threads), '-o', f'{sample}_trna_out.txt', f'{fasta}'])

    # Read the tRNAscan-SE output into a DataFrame
    df = pd.read_csv('{sample}_trna_out.txt', sep='\\t', skiprows=2)

    # Process the DataFrame
    processed_df = df[['Name', 'tRNA #', 'Begin', 'End', 'Type', 'Codon', 'Score']].drop_duplicates(subset=['Begin'])

    # Save the processed DataFrame to a new TSV file
    processed_df.to_csv('{sample}_processed_trnas.tsv', sep='\\t', index=False)
    EOF
    """


}
