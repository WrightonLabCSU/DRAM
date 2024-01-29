process COMBINE_DISTILL {

    errorStrategy 'finish'

    input:
    // No inputs for now

    output:
    path("combined.txt") into ch_combined_distill_out, optional: true

    script:
    """
    combinedChannel=""

    if [ "!{params.distill_topic}" != "" ]; then
        distill_carbon=0
        distill_energy=0
        distill_misc=0
        distill_nitrogen=0
        distill_transport=0

        validTopics=('default' 'carbon' 'energy' 'misc' 'nitrogen' 'transport')
        topics=(!{params.distill_topic.split()})

        for topic in "${topics[@]}"; do
            if [[ ! " ${validTopics[@]} " =~ " $topic " ]]; then
                echo "Invalid distill topic: $topic. Valid values are ${validTopics[@]}"
                exit 1
            fi

            case "$topic" in
                "default")
                    distill_carbon=1
                    distill_energy=1
                    distill_misc=1
                    distill_nitrogen=1
                    distill_transport=1
                    ;;
                "carbon")
                    distill_carbon=1
                    ;;
                "energy")
                    distill_energy=1
                    ;;
                "misc")
                    distill_misc=1
                    ;;
                "nitrogen")
                    distill_nitrogen=1
                    ;;
                "transport")
                    distill_transport=1
                    ;;
            esac
        done

        if [ "$distill_carbon" -eq 1 ]; then
            carbonFile="!{params.distill_carbon_sheet}"
            if [ -e "$carbonFile" ]; then
                combinedChannel+="$carbonFile,"
            else
                echo "Error: If using --distill_topic carbon (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets."
                exit 1
            fi
        fi

        if [ "$distill_energy" -eq 1 ]; then
            energyFile="!{params.distill_energy_sheet}"
            if [ -e "$energyFile" ]; then
                combinedChannel+="$energyFile,"
            else
                echo "Error: If using --distill_topic energy (or 'default'), you must have the preformatted distill sheets in ./assets/forms/distill_sheets."
                exit 1
            fi
        fi

        # Add similar conditions for distill_misc, distill_nitrogen, distill_transport...
    fi

    if [ "!{params.distill_ecosystem}" != "" ]; then
        distill_eng_sys=0
        distill_ag=0

        ecosystemList=(!{params.distill_ecosystem.split()})

        for ecosysItem in "${ecosystemList[@]}"; do
            if [[ ! " eng_sys ag " =~ " $ecosysItem " ]]; then
                echo "Invalid distill ecosystem: $ecosysItem. Valid values are eng_sys, ag"
                exit 1
            fi

            case "$ecosysItem" in
                "ag")
                    distill_ag=1
                    ;;
                "eng_sys")
                    distill_eng_sys=1
                    ;;
            esac
        done

        if [ "$distill_eng_sys" -eq 1 ]; then
            engSysFile="!{params.distill_eng_sys_sheet}"
            if [ -e "$engSysFile" ]; then
                combinedChannel+="$engSysFile,"
            else
                echo "Error: If using --distill_ecosystem eng_sys, you must have the preformatted distill sheets in ./assets/forms/distill_sheets."
                exit 1
            fi
        fi

        if [ "$distill_ag" -eq 1 ]; then
            agFile="!{params.distill_ag_sheet}"
            if [ -e "$agFile" ]; then
                combinedChannel+="$agFile,"
            else
                echo "Error: If using --distill_ecosystem ag, you must have the preformatted distill sheets in ./assets/forms/distill_sheets."
                exit 1
            fi
        fi
    fi

    if [ "!{params.distill_custom}" != "" ]; then
        customFiles=(!{params.distill_custom.replaceAll(/"/, '').split()})

        for customFile in "${customFiles[@]}"; do
            fileObject="!{customFile}"
            if [ -e "$fileObject" ]; then
                combinedChannel+="$fileObject,"
            else
                echo "Error: If using --distill_custom $customFile, you must provide the file. The path $customFile is not valid."
                exit 1
            fi
        done
    fi

    # Remove the trailing comma
    combinedChannel=${combinedChannel%,}

    echo ${combinedChannel} > combined.txt
    """
}
