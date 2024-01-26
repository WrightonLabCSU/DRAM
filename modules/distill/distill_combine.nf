process COMBINE_DISTILL {

    input:
    // Define input channels
    val ch_distill_topic
    val ch_distill_ecosys
    val ch_distill_custom

    output:
    // Define output channel
    val ch_combined_distill, emit ch_combined_distill

    script:
    // Initialize combined channel
    def combinedChannel = Channel.empty()

    // Check and add to combined channel if ch_distill_topic is not "empty"
    if (ch_distill_topic != "empty") {
        combinedChannel = combinedChannel.combine(ch_distill_topic)
    }

    // Check and add to combined channel if ch_distill_ecosys is not "empty"
    if (ch_distill_ecosys != "empty") {
        combinedChannel = combinedChannel.combine(ch_distill_ecosys)
    }

    // Check and add to combined channel if ch_distill_custom is not "empty"
    if (ch_distill_custom != "empty") {
        combinedChannel = combinedChannel.combine(ch_distill_custom)
    }

    // Set the combined channel as the output
    ch_combined_distill = combinedChannel
}
