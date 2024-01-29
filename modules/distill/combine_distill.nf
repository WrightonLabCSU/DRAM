process COMBINE_DISTILL {

    errorStrategy 'finish'

    input:
    // No inputs for now

    output:
    path("combined.txt"), emit: ch_combined_distill_out, optional: true

    script:
    """
    #!/usr/bin/env python

    import os
    import shlex

    combinedChannel = []

    print("Distill Topic:", "${params.distill_topic}")
    print("Distill Ecosystem:", "${params.distill_ecosystem}")
    print("Distill Custom:", "${params.distill_custom}")

    if "${params.distill_topic}" != "":
        distill_carbon = 0
        distill_energy = 0
        distill_misc = 0
        distill_nitrogen = 0
        distill_transport = 0

        validTopics = ['default', 'carbon', 'energy', 'misc', 'nitrogen', 'transport']
        topics = "${params.distill_topic}".split()

        print("Topics:", topics)

        for topic in topics:
            if topic not in validTopics:
                print(f"Invalid distill topic: {topic}. Valid values are {', '.join(validTopics)}")
                exit(1)

            print("Processing topic:", topic)

            if topic == 'default':
                distill_carbon = 1
                distill_energy = 1
                distill_misc = 1
                distill_nitrogen = 1
                distill_transport = 1
            elif topic == 'carbon':
                distill_carbon = 1
            elif topic == 'energy':
                distill_energy = 1
            elif topic == 'misc':
                distill_misc = 1
            elif topic == 'nitrogen':
                distill_nitrogen = 1
            elif topic == 'transport':
                distill_transport = 1

        print("Distill Flags - Carbon:", distill_carbon, "Energy:", distill_energy, "Misc:", distill_misc, "Nitrogen:", distill_nitrogen, "Transport:", distill_transport)

        if distill_carbon == 1:
            carbonFile = "${params.distill_carbon_sheet}"
            print("Carbon File:", carbonFile)
            if os.path.exists(carbonFile):
                combinedChannel.append(carbonFile)
            else:
                print(f"Error: If using --distill_topic carbon (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)

        if distill_energy == 1:
            energyFile = "${params.distill_energy_sheet}"
            if os.path.exists(energyFile):
                combinedChannel.append(energyFile)
            else:
                print(f"Error: If using --distill_topic energy (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)
        if distill_misc == 1:
            miscFile = "${params.distill_misc_sheet}"
            if os.path.exists(miscFile):
                combinedChannel.append(miscFile)
            else:
                print(f"Error: If using --distill_topic misc (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)
        if distill_nitrogen == 1:
            nitrogenFile = "${params.distill_nitrogen_sheet}"
            if os.path.exists(nitrogenFile):
                combinedChannel.append(nitrogenFile)
            else:
                print(f"Error: If using --distill_topic nitrogen (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)
        if distill_transport == 1:
            transportFile = "${params.distill_transport_sheet}"
            if os.path.exists(transportFile):
                combinedChannel.append(transportFile)
            else:
                print(f"Error: If using --distill_topic transport (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)

    if "${params.distill_ecosystem}" != "":
        distill_eng_sys = 0
        distill_ag = 0

        ecosystemList = "${params.distill_ecosystem}".split()

        for ecosysItem in ecosystemList:
            if ecosysItem not in ['eng_sys', 'ag']:
                print(f"Invalid distill ecosystem: {ecosysItem}. Valid values are eng_sys, ag")
                exit(1)

            if ecosysItem == 'ag':
                distill_ag = 1
            elif ecosysItem == 'eng_sys':
                distill_eng_sys = 1

        if distill_eng_sys == 1:
            engSysFile = "${params.distill_eng_sys_sheet}"
            if os.path.exists(engSysFile):
                combinedChannel.append(engSysFile)
            else:
                print(f"Error: If using --distill_ecosystem eng_sys, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)

        if distill_ag == 1:
            agFile = "${params.distill_ag_sheet}"
            if os.path.exists(agFile):
                combinedChannel.append(agFile)
            else:
                print(f"Error: If using --distill_ecosystem ag, you must have the preformatted distill sheets in ./assets/forms/distill_sheets.")
                exit(1)

    if "${params.distill_custom}" != "":
        # Use shlex to properly split paths considering quotes
        customFiles = shlex.split("${params.distill_custom}")

        for customFile in customFiles:
            # Ensure customFile is declared outside of the loop
            fileObject = f"{customFile}"
            if os.path.exists(fileObject):
                combinedChannel.append(fileObject)
            else:
                print(f"Error: If using --distill_custom {customFile}, you must provide the file. The path {customFile} is not valid.")
                exit(1)

    # Combine all channels into a single string
    combinedChannelStr = ','.join(combinedChannel)

    # Remove the trailing comma
    combinedChannelStr = combinedChannelStr.rstrip(',')

    # Save the combined string to combined.txt
    with open("combined.txt", "w") as output_file:
        output_file.write(combinedChannelStr)

    """
}
