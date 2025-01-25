import click

# TODO move to cloup for better constraint management
import cloup
from cloup import constraint, option, option_group
from cloup.constraints import (
    mutually_exclusive,
    If,
    RequireExactly,
    RequireAtLeast,
    AcceptAtMost,
    AnySet,
    IsSet,
    accept_none,
)

import sys
import shutil
import subprocess
import datetime as dt
from pathlib import Path

from dataclasses import dataclass

from bin.definitions import DRAM_DIR
from bin.config import Config

# from format_kegg_database import main as format_kegg_db
# from combine_annotations import main as combine_annotations
# from count_annotations import main as count_annotations
# from kofam_hmm_formatter import main as hmm_search_kofam
# from dbcan_hmm_formatter import main as hmm_search_dbcan
# from vog_hmm_formatter import main as hmm_search_vog
# from sulfur_hmm_formatter import main as hmm_search_sulfur
# from camper_hmm_formatter import main as hmm_search_camper
# from fegenie_hmm_formatter import main as hmm_search_fegenie
# from canthyd_hmm_formatter import main as hmm_search_canthyd
# from sql_add_descriptions import main as sql_add_descriptions
# from parse_annotations import main as parse_annotations
# from generate_gff_genbank import main as generate_gff_genbank
# from database
from bin.utils import split_shell_command
import bin.format_kegg_database as format_kegg_database
import bin.rename
import bin.call


def main(**kwargs):
    """
    ===========================================
                       DRAM
    (==(     )==)                 (==(     )==)
     `-.`. ,',-'                   `-.`. ,',-'
        _,-'"                         _,-'"
     ,-',' `.`-.                    ,-',' `.`-.
    (==(     )==)                  (==(     )==)
     `-.`. ,',-'                    `-.`. ,',-'
        _,-'"                         _,-'"
     ,-',' `.`-.                   ,-',' `.`-.
    (==(     )==)                 (==(     )==)
    ===========================================

    MIT License. Micobial Ecosystems Lab, Colorado State University Fort Collins. 2024 (last updated 2024)

    Description:
        The purpose of DRAM is to provide FASTA annotation, across a vast array of databases, with expertly-currated distillation.
        DRAM can be used to call, annotate and distill annotations from input FASTA files.
        Call, annotate and distill can be run together or, each can be run idependently.
    """

    config = Config(**kwargs)
    if config.format_kegg:
        format_kegg(config)

    if config.rename:
        bin.rename.main(input_fastas=config.input_fasta, output_dir=config.output_dir)
        config.set_fastas_loc(config.output_dir / "processed_inputs/*.f*")


def call(config: Config):
    bin.call.main(
        fastas=config.input_fasta,
        output_dir=config.output_dir,
        min_contig_len=config.min_contig_len,
        prodigal_mode=config.prodigal_mode,
        prodigal_trans_table=config.prodigal_trans_table
    )
    


def format_kegg(config: Config):
    """Format KEGG database for use in DRAM. Standalone operation, will exit after completion."""
    click.echo("Formatting KEGG database")
    if config.kegg_pep_root_dir:

        p = subprocess.run(
            [
                "bash",
                str(DRAM_DIR / "database/cat_kegg_pep.sh"),
                config.kegg_pep_root_dir,
                config.output_dir,
            ],
            capture_output=True,
        )
        for line in p.stdout.decode().splitlines():
            click.echo(line)
        for line in p.stderr.decode().splitlines():
            click.echo(line, err=True)
        if p.returncode != 0:
            raise click.exceptions.ClickException("Error formatting KEGG database.")
        config.kegg_pep_loc = Path(config.output_dir) / "kegg.pep"

    format_kegg_database.prepare_databases(
        kegg_loc=config.kegg_pep_loc,
        gene_ko_link_loc=config.gene_ko_link_loc,
        output_dir=config.kegg_db,
        kegg_download_date=config.kegg_download_date,
        threads=config.threads,
        skip_gene_ko_link=config.skip_gene_ko_link,
    )


if __name__ == "__main__":
    main()
