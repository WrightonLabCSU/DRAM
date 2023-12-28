process DISTILL_FINAL {

    input:
    path( metabolism_summary )

    output:
    path( "distillate.xlsx" ), emit: distillate

    script:
    """

    python /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/assets/generate_multi_sheet_xlsx.py --input-file ${metabolism_summary} --output-file distillate.xlsx


    """
}
