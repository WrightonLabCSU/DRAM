process COMBINE_DISTILL {

    input:
    // Define input channels
    path( ch_distill_carbon, stageAs:  "carbon_distill_sheet.tsv")
    path( ch_distill_energy, stageAs:  "energy_distill_sheet.tsv") 
    path( ch_distill_misc, stageAs:  "misc_distill_sheet.tsv") 
    path( ch_distill_nitrogen, stageAs:  "nitrogen_distill_sheet.tsv") 
    path( ch_distill_transport, stageAs:  "transport_distill_sheet.tsv") 
    path( ch_distill_ag, stageAs:  "ag_distill_sheet.tsv") 
    path( ch_distill_eng_sys, stageAs:  "eng_sys_distill_sheet.tsv*")
    path( ch_distill_custom, stageAs:  "custom_distill_sheet.tsv") 

    output:
    tuple path( ch_distill_carbon ), path( ch_distill_energy ), path( ch_distill_misc ), path( ch_distill_nitrogen ), path( ch_distill_transport ), path( ch_distill_ag ), path( ch_distill_eng_sys), path( ch_distill_custom ), emit: ch_combined_distill_sheets


    script:
    """
    """
}
