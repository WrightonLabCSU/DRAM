process RENAME_FASTA {


    tag { sample }

    input:
    tuple val(sample), path(fasta)
    
    output:
    tuple val(sample), path("*.fna"), emit: renamed_fasta


    script:

    """

    rename.sh \\
    in=${fasta} \\
    out=${sample}_renamed.fna \\
    prefix=${sample} \\
    addprefix=t

    """
}  