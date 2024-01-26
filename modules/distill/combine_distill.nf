process COMBINE_DISTILL {

    input:
    // Define input channels
    val ch_distill_topic
    val ch_distill_ecosys
    val ch_distill_custom

    output:
    path("combined.txt"), emit: ch_combined_distill_out

    shell:
    '''
    combinedChannel=""

    # Debug statements within COMBINE_DISTILL shell block
    echo "Debug: ch_distill_topic = !{ch_distill_topic}"
    echo "Debug: ch_distill_ecosys = !{ch_distill_ecosys}"
    echo "Debug: ch_distill_custom = !{ch_distill_custom}"
    
    # Check and add to combined channel if ch_distill_ecosys is not "empty"
    if [[ "!{ch_distill_topic}" != "empty" ]]; then
        combinedChannel="${combinedChannel}!{ch_distill_topic},"
    fi

    # Check and add to combined channel if ch_distill_topic is not "empty"
    if [[ "!{ch_distill_ecosys}" != "empty" ]]; then
        combinedChannel="${combinedChannel}!{ch_distill_ecosys},"
    fi

    # Check and add to combined channel if ch_distill_custom is not "empty"
    if [[ "!{ch_distill_custom}" != "empty" ]]; then
        combinedChannel="${combinedChannel}!{ch_distill_custom},"
    fi

    echo $combinedChannel > combined.txt

    # Print debug information about the channels after concatenation
    echo "Debug: combinedChannel = ${combinedChannel}"


    # Print the contents of the resulting combined.txt file
    cat combined.txt
    '''
}
