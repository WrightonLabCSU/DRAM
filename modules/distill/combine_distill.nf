process COMBINE_DISTILL {

    input:
    // Define input channels

    val ch_combine_test

    //val ch_distill_topic
    //val ch_distill_ecosys
    //val ch_distill_custom
    //val distill_flag_real

    output:
    //path("combined.txt"), emit: ch_combined_distill_out

    //when:
    //distill_flag_real == "1"

    script:
    """
    echo ${ch_combine_test}
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