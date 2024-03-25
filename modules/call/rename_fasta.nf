process RENAME_FASTA {
    conda = './assets/conda/environment.yml'

    tag { sample }

    input:
    tuple val(sample), path(fasta)
    
    output:
    tuple val(sample), path("*.fna"), emit: renamed_fasta


    script:

    """
    conda activate /home/rwoyda/miniconda3/envs/dram2-env
    echo "Active Python version:"
    python --version
    echo "Which Python:"
    which python
    which conda
    echo "List conda envs"
    conda env list
    echo "Which conda"
    which conda
    rename.sh \\
    in=${fasta} \\
    out=${sample}_renamed.fna \\
    prefix=${sample} \\
    addprefix=t

    """
}  