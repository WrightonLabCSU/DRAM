process RENAME_FASTA {

    tag { "renaming_fastas" }

    conda "${moduleDir}/environment.yml"

    input:
    val fasta_names
    path fastas

    output:
    path("renamed/*.fna"), emit: renamed_fasta_paths

    script:
    // def samples = fasta_names as List

    def rename_cmds = fasta_names.indices.collect { i ->
        def name = fasta_names[i]
        def file = fastas[i]
        "rename.sh in=${file} out=renamed/${name}.fna prefix=${name} addprefix=t"
    }.join("\n")


    """
    mkdir -p renamed
    ${rename_cmds}
    """
}
