process COMBINE_DISTILL {

    input:
    val queueChannel

    output:
    path("combined.txt"), emit: ch_combined_distill_out

    script:
    """
    # Dequeue the information from the queue channel
    info=($queueChannel)
    ch_distill_topic=${info[0]}
    ch_distill_ecosys=${info[1]}
    ch_distill_custom=${info[2]}
    
    # Your existing script for combining the channels
    combinedChannel=""
    if [[ "${ch_distill_topic}" != "empty" ]]; then
        combinedChannel="${combinedChannel}${ch_distill_topic},"
    fi
    if [[ "${ch_distill_ecosys}" != "empty" ]]; then
        combinedChannel="${combinedChannel}${ch_distill_ecosys},"
    fi
    if [[ "${ch_distill_custom}" != "empty" ]]; then
        combinedChannel="${combinedChannel}${ch_distill_custom},"
    fi
    echo $combinedChannel > combined.txt
    """
}