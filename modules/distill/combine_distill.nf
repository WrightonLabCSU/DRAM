process COMBINE_DISTILL {

    input:
    // Define input channels
    val ch_distill_topic
    val ch_distill_ecosys
    val ch_distill_custom

    output:
    // Define output channel

    script:
    """
    echo "ch_distill_topic: ${ch_distill_topic}"
    echo "ch_distill_ecosys: ${ch_distill_ecosys}"
    echo "ch_distill_custom: ${ch_distill_custom}"

    combinedChannel=""

    # Check and add to combined channel if ch_distill_topic is not "empty"
    if [[ "${ch_distill_ecosys}" != "empty" ]]; then
        combinedChannel="${combinedChannel}${ch_distill_ecosys}"
    fi

    # Check and add to combined channel if ch_distill_ecosys is not "empty"
    if [[ "${ch_distill_topic}" != "empty" ]]; then
        combinedChannel="${combinedChannel}${ch_distill_topic}"
    fi

    # Check and add to combined channel if ch_distill_custom is not "empty"
    if [[ "${ch_distill_custom}" != "empty" ]]; then
        combinedChannel="${combinedChannel}${ch_distill_custom}"
    fi

    echo $combinedChannel

    """
}
