from pathlib import Path
import subprocess

from bin.utils import split_shell_command

def main(input_fastas: Path, output_dir: Path) -> None:
    for fasta in ((input_fastas.parent)  # The parent directory of the input path
                  .glob(input_fastas.name)):  # Glob for all files matching the input pattern
        sample: str = fasta.stem.replace(".", "_")
        subprocess.run(split_shell_command(f"bash rename.sh in={fasta} out={output_dir}/processed_inputs/{sample}_renamed.fna prefix=${sample} addprefix=t"))
