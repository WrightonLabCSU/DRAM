include { logColours           } from '../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion           } from '../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Print pipeline summary on completion
//
def getFastaChannel(input_fasta, fasta_fmt) {
    
    def ch_fasta = Channel
        .fromPath(file(input_fasta) / fasta_fmt, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find any fasta files matching: ${input_fasta}\nNB: Path needs to follow pattern: path/to/directory/" }
    
    ch_fasta = ch_fasta.map {
        sampleName = it.getName().replaceAll(/\.[^.]+$/, '').replaceAll(/\./, '-')
        tuple(sampleName, it)
    }
    return ch_fasta  // channel: [ val(sample name), path(fasta) ]
}

//
// DRAM logo
//
def dramLogo(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    String.format(  // In groovy we can use $/ to start a string to avoid escaping (end with /$)
        $/\n
        =====================================================================        
                         _____    _____               __  __            
        (==(     )==)   |  __ \  |  __ \      /\     |  \/  |   (==(     )==)    
         `-.`. ,',-'    | |  | | | |__) |    /  \    | \  / |    `-.`. ,',-'     
            _,-'"       | |  | | |  _  /    / /\ \   | |\/| |       _,-'"        
         ,-',' `.`-.    | |__| | | | \ \   / ____ \  | |  | |    ,-',' `.`-.     
        (==(     )==)   |_____/  |_|  \_\ /_/    \_\ |_|  |_|   (==(     )==)    
                                            
        =====================================================================
        /$.stripIndent()  // end the string with /$ to avoid escaping
    )
}

//
// Citation string for pipeline  TODO: finish citation information
//
def workflowCitation() {
    def temp_doi_ref = ""
    String[] manifest_doi = workflow.manifest.doi.tokenize(",")
    // Using a loop to handle multiple DOIs
    // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
    // Removing ` ` since the manifest.doi is a string and not a proper list
    for (String doi_ref: manifest_doi) temp_doi_ref += "  https://doi.org/${doi_ref.replace('https://doi.org/', '').replace(' ', '')}\n"
    return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* DRAM\n" +
        temp_doi_ref + "\n" +
        // "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
        "* Software dependencies\n" +
        "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
}  // TODO add DOI for DRAM