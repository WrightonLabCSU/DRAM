//
// Subworkflow that uses the nf-schema plugin to validate parameters and render the parameter summary
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    take:
    input_workflow      // workflow: the workflow object used by nf-schema to get metadata from the workflow
    validate_params     // boolean:  validate the parameters
    parameters_schema   // string:   path to the parameters JSON schema.
                        //           this has to be the same as the schema given to `validation.parametersSchema`
                        //           when this input is empty it will automatically use the configured schema or
                        //           "${projectDir}/nextflow_schema.json" as default. This input should not be empty
                        //           for meta pipelines

    main:

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //
    if(parameters_schema) {
        log.info paramsSummaryLog(input_workflow, parameters_schema:parameters_schema)
    } else {
        log.info paramsSummaryLog(input_workflow)
    }

    //
    // Coerce integer params that arrive as strings from the CLI.
    // nf-schema 2.6.x lenientMode coerces booleans but not numerics.
    //
    [
        'array_size', 'queue_size', 'threads',
        'tiny_cpus_limit', 'small_cpus_limit', 'medium_cpus_limit', 'big_cpus_limit', 'huge_cpus_limit',
        'tiny_gb_mem_limit', 'small_gb_mem_limit', 'medium_gb_mem_limit', 'big_gb_mem_limit', 'huge_gb_mem_limit',
        'tiny_hr_time_limit', 'small_hr_time_limit', 'medium_hr_time_limit', 'big_hr_time_limit', 'huge_hr_time_limit',
        'cpu_provision_limit', 'mem_gb_provision_limit', 'time_hr_provision_limit',
        'kofam_chunk_size', 'vog_chunk_size', 'amg_length_from_end',
    ].each { k ->
        if (params[k] instanceof String) params[k] = params[k] as Integer
    }

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if(validate_params) {
        if(parameters_schema) {
            validateParameters(parameters_schema:parameters_schema)
        } else {
            validateParameters()
        }
    }

    emit:
    dummy_emit = true
}

