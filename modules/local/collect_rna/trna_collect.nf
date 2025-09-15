#!/usr/bin/env nextflow

process TRNA_COLLECT {
    label 'process_small'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_pandas_barrnap_trnascan-se:ed2ab26abf39304b"

    input:
    file combined_trnas

    output:
    path("collected_trnas.tsv"), emit: trna_collected_out, optional: true
    path("combined_trna_scan.tsv"), emit: trna_combined_out, optional: true

    script:
    """
    # export constants for script
    export FASTA_COLUMN="${params.CONSTANTS.FASTA_COLUMN}"

    awk '(NR==1) || (FNR>1)' *processed_trnas.tsv > combined_trna_scan.tsv

    trna_collect.py
    """
}
