from pathlib import Path
import subprocess

from bin.utils import split_shell_command, get_fasta_sample_name

def main(fastas: Path, output_dir: Path) -> None:
    for fasta in ((fastas.parent)  # The parent directory of the input path
                  .glob(fastas.name)):  # Glob for all files matching the input pattern
        sample = get_fasta_sample_name(fasta)
        subprocess.run(split_shell_command(f"bash rename.sh in={fasta} out={output_dir}/processed_inputs/{sample}_renamed.fna prefix=${sample} addprefix=t"))
