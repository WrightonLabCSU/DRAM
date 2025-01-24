from dataclasses import dataclass, field
from typing import List, Optional
from pathlib import Path
import datetime as dt


@dataclass
class Config:
    # General Options
    
    output_dir: Path = field(default_factory=lambda: Path.cwd() / f"DRAM_output.{dt.datetime.now().strftime('%Y%m%dT%H%M%S')}")
    input_fasta: Path = field(default_factory=lambda: Path.cwd() / "input_fasta")
    fasta_fmt: str = "*.f*"
    fastas: Path = None
    rename: bool = False
    threads: int = 10

    # Call Options
    call: Optional[bool] = None
    prodigal_mode: str = "single"
    prodigal_trans_table: int = 1
    min_contig_len: int = 2500

    # Annotate Options
    annotate: Optional[bool] = None
    use_db: List[str] = field(default_factory=list)
    use_dbset: Optional[str] = None
    input_genes: Optional[Path] = None
    add_annotations: Optional[Path] = None
    generate_gff: bool = False
    generate_gbk: bool = False

    # Distill Options
    distill: Optional[bool] = None
    annotations: Optional[Path] = None
    rrnas: Optional[Path] = None
    trnas: Optional[Path] = None
    bin_quality: Optional[Path] = None
    taxa: Optional[Path] = None
    distill_topic: Optional[str] = None
    distill_ecosystem: Optional[str] = None
    distill_custom: Optional[Path] = None

    # Product Options
    product: Optional[bool] = None
    groupby_column: str = "fasta"
    module_steps_form: Optional[Path] = None
    etc_module_form: Optional[Path] = None
    function_heatmap_form: Optional[Path] = None

    # Format KEGG Options
    format_kegg: Optional[bool] = None
    kegg_pep_loc: Optional[Path] = None
    kegg_pep_root_dir: Optional[Path] = None
    gene_ko_link_loc: Optional[Path] = None
    kegg_db: Path = field(default_factory=lambda: Path.cwd() / "kegg_db")
    kegg_download_date: Optional[str] = None
    skip_gene_ko_link: bool = False


    def __post_init__(self):
        if not self.fastas:  # if fastas is not set, infer from input_fasta and fasta_fmt
            self.fastas = Path(self.input_fasta) / self.fasta_fmt
        # ensure input_fasta and fasta_fmt are derived from fastas and not the defaults
        self.input_fasta = self.fastas.parent
        self.fasta_fmt = self.fastas.name