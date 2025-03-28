process RENAME_FASTA {


    tag { input_fasta }

    input:
    tuple val(input_fasta), path(fasta)
    
    output:
    tuple val(input_fasta), path("*.fna"), emit: renamed_fasta


    script:

    """

    rename.sh \\
    in=${fasta} \\
    out=${input_fasta}_renamed.fna \\
    prefix=${input_fasta} \\
    addprefix=t

    """
}  