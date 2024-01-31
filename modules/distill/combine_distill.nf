process COMBINE_DISTILL {

    input:
    // Define input channels
    path( ch_distill_carbon, stageAs:  "carbon.tsv")
    path( ch_distill_energy, stageAs:  "energy.tsv") 
    path( ch_distill_misc, stageAs:  "misc.tsv") 
    path( ch_distill_nitrogen, stageAs:  "nitrogen.tsv") 
    path( ch_distill_transport, stageAs:  "transport.tsv") 
    path( ch_distill_ag, stageAs:  "ag.tsv") 
    path( ch_distill_eng_sys, stageAs:  "eng_sys.tsv*")
    path( ch_distill_custom, stageAs:  "custom.tsv") 

    output:
    tuple path( ch_distill_carbon ), path( ch_distill_energy ), path( ch_distill_misc ), path( ch_distill_nitrogen ), path( ch_distill_transport ), path( ch_distill_ag ), path( ch_distill_eng_sys), path( ch_distill_custom ), emit: ch_combined_distill_sheets


    script:
    """


    """
}

/*
    shell:
    '''
    combinedChannel=""
    
    # Check and add to combined channel if ch_distill_topic is not "empty"
    if [[ "!{ch_distill_ecosys}" != "empty" ]]; then
        combinedChannel="${combinedChannel}!{ch_distill_ecosys},"
    fi

    # Check and add to combined channel if ch_distill_ecosys is not "empty"
    if [[ "!{ch_distill_topic}" != "empty" ]]; then
        combinedChannel="${combinedChannel}!{ch_distill_topic},"
    fi

    # Check and add to combined channel if ch_distill_custom is not "empty"
    if [[ "!{ch_distill_custom}" != "empty" ]]; then
        combinedChannel="${combinedChannel}!{ch_distill_custom},"
    fi

    echo $combinedChannel > combined.txt
    '''

    */