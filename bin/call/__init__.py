from pathlib import Path
from os import stat
import logging

from bin.utils import run_process, get_fasta_sample_name

logger = logging.getLogger(__name__)

def main(fastas: Path, output_dir: Path, min_contig_len: int, prodigal_mode: str, prodigal_trans_table: int) -> None:
    for fasta in ((fastas.parent)  # The parent directory of the input path
                  .glob(fastas.name)):  # Glob for all files matching the input pattern
        call(output_dir, min_contig_len, fasta, prodigal_mode, prodigal_trans_table)

def call(output_dir, min_contig_len, fasta, prodigal_mode, prodigal_trans_table):
    sample = get_fasta_sample_name(fasta)
    prodigal_dir = Path(output_dir) / "Prodigal_v2.6.3"
    run_process(
            f"bash reformat.sh in={fasta} out={prodigal_dir / f'{sample}_{min_contig_len}.fa'} minlength={min_contig_len}"
        )
    if stat(prodigal_dir / f'{sample}_{min_contig_len}.fa').st_size == 0:
        logger.warning(f"Sample {sample} - Error: no contigs after filtering for minimum contig length ({min_contig_len})")
        return
    
    run_process(
        f"prodigal -i {prodigal_dir / sample}_{min_contig_len}.fa -o {prodigal_dir / sample}_called_genes.gff -p {prodigal_mode} -g {prodigal_trans_table} -f gff -d {prodigal_dir / sample}_called_genes.fna -a {prodigal_dir / sample}_called_genes.faa"
    )
    if stat(f"{prodigal_dir / sample}_called_genes.gff").st_size == 0:
        logger.warning(f"Sample {sample} - Warning: called genes GFF file is empty or does not exist.")
    