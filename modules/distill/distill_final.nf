process DISTILL_FINAL {

    input:
    path( metabolism_summary )
    file( ch_distill_final_script )

    output:
    path( "distillate.xlsx" ), emit: distillate

    script:
    """

    python ${ch_distill_final_script} --input-file ${metabolism_summary} --output-file distillate.xlsx

    """
}
