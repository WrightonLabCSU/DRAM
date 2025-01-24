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

from pathlib import Path

from bin.main import main as dram


@cloup.command(help=dram.__doc__, show_constraints=True)
@cloup.option_group(
    "General Options",
    cloup.option(
        "--output_dir",
        type=click.Path(path_type=Path),
        default="./DRAM_output.{datetime}",
        help="Output directory",
        show_default=True,
    ),
    cloup.option(
        "--input_fasta",
        type=click.Path(path_type=Path),
        default="./input_fasta",
        help="Directory containing input FASTA files. Used in multiple pipeline steps",
        show_default=True,
    ),
    cloup.option(
        "--fasta_fmt",
        type=click.STRING,
        default="*.f*",
        help="Glob format for matching fastas",
        show_default=True,
    ),
    cloup.option(
        "--rename",
        is_flag=True,
        help=(
            "Rename FASTA headers based on file name."
            " Example: sample1.fa --> (header renamed to) >sample1. "
            " Why? DRAM output is focused on scaffolds/contigs with respect to each provided input sample."
            " Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguashable."
            " *If you have already renamed your FASTA headers, do not include this flag."
        ),
    ),
    cloup.option(
        "--threads",
        type=click.INT,
        default=10,
        help="Number of threads to use for processing.",
        show_default=True,
    ),
)
@cloup.option_group(
    "Call Options",
    "Call Description: The purpose of DRAM call is to call genes on input FASTA files using Prodigal.",
    cloup.option(
        "--call", help="Run `Call` genes on the input FASTA files using Prodigal."
    ),
    cloup.option(
        "--prodigal_mode",
        type=click.Choice(["single", "meta"]),
        default="single",
        help="Prodigal mode for gene calling. Default: 'single'.",
    ),
    cloup.option(
        "--prodigal_trans_table",
        type=click.IntRange(1, 25),
        default=1,
        help="Translation table to use (1-25). Default: '1'.",
    ),
    cloup.option(
        "--min_contig_len",
        type=click.INT,
        default=2500,
        help="Minimum contig length in base pairs. Default: '2500'.",
    ),
)
@cloup.option_group(
    "Annotate Options",
    "Annotate Description: The purpose of DRAM '--annotate' is to annotate called genes on input (nucleotide) FASTA (fa*) files using downloaded databases.",
    cloup.option(
        "--annotate", help="Run `Annotate` called genes using downloaded databases."
    ),
    cloup.option(
        "--use_db",
        multiple=True,
        type=click.Choice(
            [
                "camper",
                "cant_hyd",
                "dbcan",
                "fegenie",
                "kegg",
                "kofam",
                "merops",
                "methyl",
                "heme",
                "pfam",
                "sulfur",
                "uniref",
            ]
        ),
        help="Specify databases to use. Multiple values allowed.",
    ),
    cloup.option(
        "--use_dbset",
        type=click.Choice(
            [
                "metabolism_kegg_set",
                "metabolism_set",
                "adjectives_kegg_set",
                "adjectives_set",
            ]
        ),
        help=(
            "Use a predefined database set. "
            "Choices: metabolism_kegg_set, metabolism_set, adjectives_kegg_set, adjectives_set."
        ),
    ),
    cloup.option(
        "--input_genes",
        type=click.Path(path_type=Path),
        help="Directory containing called genes (.faa).",
    ),
    cloup.option(
        "--add_annotations",
        type=click.Path(path_type=Path),
        help="Path to old annotations to include in the current run.",
    ),
    cloup.option(
        "--generate_gff",
        is_flag=True,
        help="Generate output GFF files for each sample.",
    ),
    cloup.option(
        "--generate_gbk",
        is_flag=True,
        help="Generate output GBK files for each sample.",
    ),
)
# The following will be shown (with --help) under "Other options"
@cloup.option_group(
    "Distill Options",
    (
        "Distill Description: The purpose of DRAM --distill is to distill down annotations based on curated distillation summary form(s)."
        " User's may also provide a custom distillate via --distill_custom <path/to/file> (TSV forms)."
        " Distill can be ran independent of --call and --annotate however, annotations must be provided (--annotations <path/to/annotations.tsv>)."
        " Optional tRNA, rRNA and bin quality may also be provided."
    ),
    cloup.option(
        "--distill",
        help="Run `Distill` to distill the annotations into a multi-sheet distillate.xlsx.",
    ),
    cloup.option(
        "--annotations",
        type=click.Path(path_type=Path),
        help="Path to annotations.tsv file. Required if not running --call and --annotate.",
    ),
    cloup.option(
        "--rrnas", type=click.Path(path_type=Path), help="Path to rRNA.tsv file."
    ),
    cloup.option(
        "--trnas", type=click.Path(path_type=Path), help="Path to tRNA.tsv file."
    ),
    cloup.option(
        "--bin_quality",
        type=click.Path(path_type=Path),
        help="Path to bin-quality.tsv file.",
    ),
    cloup.option(
        "--taxa", type=click.Path(path_type=Path), help="Path to bin-taxonomy.tsv file."
    ),
    cloup.option(
        "--distill_topic",
        type=str,
        help=(
            "Topics for distillation (e.g., 'carbon transport')."
            " If multiple topics, enclose in quotes."
        ),
    ),
    cloup.option(
        "--distill_ecosystem",
        type=str,
        help=(
            "Ecosystem for distillation (e.g., 'eng_sys ag')."
            " If multiple ecosystems, enclose in quotes."
        ),
    ),
    cloup.option(
        "--distill_custom",
        type=click.Path(path_type=Path),
        help="Path to custom distillate TSV file.",
    ),
)
@cloup.option_group(
    "Product Options",
    (
        "Product Description: The purpose of DRAM --product is to generate a product visualization of the annotations"
        " and save the output to the output directory."
    ),
    cloup.option(
        "--product",
        help="Run `Product` to generate a product visualization of the annotations and save the output to the output directory.",
    ),
    cloup.option(
        "--annotations",
        type=click.Path(path_type=Path),
        required=True,
        help="Path to annotations.tsv file.",
    ),
    cloup.option(
        "--groupby_column",
        type=str,
        default="fasta",
        help="Column for grouping in annotations file. Default: 'fasta'.",
    ),
    cloup.option(
        "--module_steps_form",
        type=click.Path(path_type=Path),
        help="Path to module steps form TSV file.",
    ),
    cloup.option(
        "--etc_module_form",
        type=click.Path(path_type=Path),
        help="Path to etc module form TSV file.",
    ),
    cloup.option(
        "--function_heatmap_form",
        type=click.Path(path_type=Path),
        help="Path to function heatmap form TSV file.",
    ),
)
@cloup.option_group(
    "Format KEGG Options",
    (
        "Format KEGG DB Description: The purpose of DRAM '--format_kegg' is to format the raw KEGG database (pep files)"
        " for use in DRAM. Because use of the KEGG database requires a license, the user must download the raw KEGG database"
        " themselves. Currently this only supports formatting a concatenated version of the KEGG pep files see"
        " (https://github.com/WrightonLabCSU/DRAM/issues/305) for guidance on how to prepare your pep files for formatting."
    ),
    cloup.option(
        "--format_kegg",
        help="Run `Format KEGG` to Format the KEGG database for use in DRAM. Standalone operation, will exit after completion.",
    ),
    cloup.option(
        "--kegg_pep_loc",
        type=click.Path(path_type=Path),
        help="Path to KEGG peptide location."
        " Either pass the kegg_pep_loc directly, or the path to the downloaded pep files from your"
        " KEGG database with `--kegg_pep_root_dir`, but not both.",
    ),
    cloup.option(
        "--kegg_pep_root_dir",
        type=click.Path(path_type=Path),
        help="Root directory path to KEGG pep files, from the downloaded pep files from your KEGG database."
        " If included, will concatenate all pep files together to one pep file. "
        " Should be in the directory format: `{kegg_pep_root_dir}/*/*.pep`",
    ),
    cloup.option(
        "--gene_ko_link_loc",
        type=click.Path(path_type=Path),
        help="Path to genes_ko_link file.",
    ),
    cloup.option(
        "--kegg_db",
        type=click.Path(path_type=Path),
        default="./kegg_db",
        help="Path to save formatted KEGG database. Default: './kegg_db'.",
    ),
    cloup.option(
        "--kegg_download_date",
        type=str,
        help="Date KEGG database was downloaded (YYYY-MM-DD).",
    ),
    cloup.option(
        "--skip_gene_ko_link", is_flag=True, help="Skip the gene_ko_link file check."
    ),
)
@constraint(
    If("format_kegg", then=AcceptAtMost(0)), ["call", "annotate", "distill", "product"]
)
@constraint(
    If("format_kegg", then=RequireExactly(1)), ["kegg_pep_loc", "kegg_pep_root_dir"]
)
def cli(**kwargs):
    dram(**kwargs)


if __name__ == "__main__":
    cli()
