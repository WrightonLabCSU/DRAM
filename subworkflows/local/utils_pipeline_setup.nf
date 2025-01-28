
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