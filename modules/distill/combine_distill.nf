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

    combined="banana"
    echo $combined



    """
}
