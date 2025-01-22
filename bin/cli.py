
import click
# TODO move to cloup for better constraint management
from cloup import constraint
from cloup.constraints import mutually_exclusive, If, RequireExactly

import sys
import shutil
import subprocess
import datetime as dt
from pathlib import Path

from bin.definitions import DRAM_DIR
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


PIPELINE_ORDER = [
    'call', 'annotate', 'distill', 'product'
]
ALL_PIPELINE_COMMANDS = ["format-kegg"] + PIPELINE_ORDER

@click.group(chain=True,
             invoke_without_command=True
             )
@click.option('--output_dir', type=click.Path(path_type=Path), default='./DRAM_output.{datetime}',
              help="Output directory", show_default=True)
@click.option('--input_fasta', type=click.Path(path_type=Path), default='./input_fasta',
              help="Directory containing input FASTA files. Used in multiple pipeline steps", show_default=True)
@click.option('--fasta_fmt', type=click.STRING, default='*.f*',
              help="Glob format for matching fastas", show_default=True)
@click.option('--rename', is_flag=True, help=(
        "Rename FASTA headers based on file name."
        " Example: sample1.fa --> (header renamed to) >sample1. "
        " Why? DRAM output is focused on scaffolds/contigs with respect to each provided input sample."
        " Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguashable."
        " *If you have already renamed your FASTA headers, do not include this flag."
))
@click.option('--threads', type=click.INT, default=10,
              help="Number of threads to use for processing.", show_default=True)
@click.pass_context
def cli(ctx, output_dir, input_fasta, fasta_fmt, rename, threads):
    """Main CLI tool."""
    # click.echo(ctx.invoked_subcommand)
    if ctx.invoked_subcommand is None and not rename:
        click.echo(ctx.get_help())
        return
    ctx.ensure_object(dict)
    
    
    if output_dir == Path('./DRAM_output.{datetime}'):
        output_dir = Path.cwd() / f'DRAM_output.{dt.datetime.now().strftime("%Y%m%dT%H%M%S")}'
        ctx.params['output_dir'] = output_dir
    
    if output_dir.exists():
        response = input("Output directory already exists. Do you want to overwrite? (y/n)")
        if response.lower() != 'y':
            sys.exit("Exiting.")
        else:
            click.echo(f"Removing existing directory: {output_dir} and creating a new one.")
            shutil.rmtree(output_dir)

    (output_dir / "processed_inputs").mkdir(parents=True, exist_ok=False)
    
    commands = sys.argv[1:]  # We can use .remove() later because slicing creates a copy

    if "format-kegg" in commands and any(cmd in PIPELINE_ORDER for cmd in commands):
        raise click.exceptions.BadArgumentUsage("format-kegg must be run alone.")
    if "format-kegg" in commands:  # if just format-kegg, we run it
        return
    
    # Else we check the rest of the steps

    pipeline_commands_in_order = [cmd for cmd in sys.argv[1:] if cmd in ALL_PIPELINE_COMMANDS]
    pipeline_commands_in_order_working = pipeline_commands_in_order.copy()

    click.echo(f"Running DRAM pipeline with commands: {pipeline_commands_in_order}")

    # First we check that there aren't duplicate pipeline commands in the CLI
    # we use sys.argv here instead of commands since we are removing later in the loop,
    # and can't remove from the loop we are iterating from
    # We use set, to make sure we hit each cmd only once, so if it is cmd1 cmd2 cmd1, we only remove cmd1 once,
    # and detect the duplicate
    for cmd in set(pipeline_commands_in_order):
        if cmd in PIPELINE_ORDER:
            pipeline_commands_in_order_working.remove(cmd)

    if pipeline_commands_in_order_working:
        raise click.exceptions.BadArgumentUsage(f"1 or more pipeline commands entered more than once in DRAM launch command: {pipeline_commands_in_order_working}")

    # Next we check to make sure the CLI pipeline commands are in the right order, not every
    # command has to be in `PIPELINE_ORDER`, the cmds that are, have to be a subset in that order
    for i in range(len(pipeline_commands_in_order) - 1):
        if PIPELINE_ORDER.index(pipeline_commands_in_order[i]) > PIPELINE_ORDER.index(pipeline_commands_in_order[i + 1]):
            raise click.exceptions.BadArgumentUsage(
                f"Pipeline commands are out of order: {pipeline_commands_in_order}. Please follow a subset of the order: {PIPELINE_ORDER}"
            )
    input_fasta = input_fasta / fasta_fmt
    ctx.params['input_fasta'] = input_fasta
    if rename:
        bin.rename.main(input_fastas=input_fasta, output_dir=output_dir)


@cli.command()
@click.option('--prodigal_mode', type=click.Choice(['single', 'meta']), default='single',
              help="Prodigal mode for gene calling. Default: 'single'.")
@click.option('--prodigal_tras_table', type=click.IntRange(1, 25), default=1,
              help="Translation table to use (1-25). Default: '1'.")
@click.option('--min_contig_len', type=int, default=2500,
              help="Minimum contig length in base pairs. Default: '2500'.")
def call(prodigal_mode, prodigal_tras_table, min_contig_len):
    """Call genes using Prodigal."""
    click.echo("Calling genes")
    # TODO: Implement corresponding Python function


@cli.command()
@click.option('--use_db', multiple=True, type=click.Choice([
    'camper', 'cant_hyd', 'dbcan', 'fegenie', 'kegg', 'kofam', 'merops', 'methyl',
    'heme', 'pfam', 'sulfur', 'uniref'
]), help="Specify databases to use. Multiple values allowed.")
@click.option('--use_dbset', type=click.Choice([
    'metabolism_kegg_set', 'metabolism_set', 'adjectives_kegg_set', 'adjectives_set'
]), help=(
        "Use a predefined database set. "
        "Choices: metabolism_kegg_set, metabolism_set, adjectives_kegg_set, adjectives_set."
))
@click.option('--input_genes', type=click.Path(path_type=Path),
              help="Directory containing called genes (.faa).")
@click.option('--add_annotations', type=click.Path(path_type=Path),
              help="Path to old annotations to include in the current run.")
@click.option('--generate_gff', is_flag=True, help="Generate output GFF files for each sample.")
@click.option('--generate_gbk', is_flag=True, help="Generate output GBK files for each sample.")
def annotate(use_db, use_dbset, input_genes, add_annotations, generate_gff, generate_gbk):
    """Annotate called genes using downloaded databases"""
    click.echo("Annotating genes")
    # combine_annotations(input_dir=input_dir, output_dir=output_dir)


@cli.command()
@click.option('--annotations', type=click.Path(path_type=Path),
              help="Path to annotations.tsv file. Required if not running --call and --annotate.")
@click.option('--rrnas', type=click.Path(path_type=Path), help="Path to rRNA.tsv file.")
@click.option('--trnas', type=click.Path(path_type=Path), help="Path to tRNA.tsv file.")
@click.option('--bin_quality', type=click.Path(path_type=Path), help="Path to bin-quality.tsv file.")
@click.option('--taxa', type=click.Path(path_type=Path), help="Path to bin-taxonomy.tsv file.")
@click.option('--distill_topic', type=str, help=(
        "Topics for distillation (e.g., 'carbon transport'). "
        "If multiple topics, enclose in quotes."))
@click.option('--distill_ecosystem', type=str, help=(
        "Ecosystem for distillation (e.g., 'eng_sys ag'). "
        "If multiple ecosystems, enclose in quotes."))
@click.option('--distill_custom', type=click.Path(path_type=Path),
              help="Path to custom distillate TSV file.")
def distill(annotations, rrnas, trnas, bin_quality, taxa, distill_topic, distill_ecosystem, distill_custom):
    """Distill the annotations into a multi-sheet distillate.xlsx"""
    click.echo("Distilling annotations")
    # count_annotations(annotations=annotations, output=output)

@cli.command()
@click.option('--annotations', type=click.Path(path_type=Path), required=True,
              help="Path to annotations.tsv file.")
@click.option('--groupby_column', type=str, default='fasta',
              help="Column for grouping in annotations file. Default: 'fasta'.")
@click.option('--module_steps_form', type=click.Path(path_type=Path),
              help="Path to module steps form TSV file.")
@click.option('--etc_module_form', type=click.Path(path_type=Path),
              help="Path to etc module form TSV file.")
@click.option('--function_heatmap_form', type=click.Path(path_type=Path),
              help="Path to function heatmap form TSV file.")
def product(annotations, groupby_column, module_steps_form, etc_module_form, function_heatmap_form):
    """Generate a product visualization of the annotations."""
    click.echo("Generating product visualization")

@cli.command()
@click.option('--kegg_pep_loc', type=click.Path(path_type=Path),
              help="Path to KEGG peptide location."
                   " Either pass the kegg_pep_loc directly, or the path to the downloaded pep files from your"
                   " KEGG database with `--kegg_pep_root_dir`, but not both.")
@click.option('--kegg_pep_root_dir', type=click.Path(path_type=Path),
              help="Root directory path to KEGG pep files, from the downloaded pep files from your KEGG database."
                   " If included, will concatenate all pep files together to one pep file. "
                   " Should be in the directory format: `{kegg_pep_root_dir}/*/*.pep`", required=False)
@click.option('--gene_ko_link_loc', type=click.Path(path_type=Path), help="Path to genes_ko_link file.")
@click.option('--kegg_db', type=click.Path(path_type=Path), default='./kegg_db',
              help="Path to save formatted KEGG database. Default: './kegg_db'.")
@click.option('--kegg_download_date', type=str, help="Date KEGG database was downloaded (YYYY-MM-DD).")
@click.option('--skip_gene_ko_link', is_flag=True, help="Skip the gene_ko_link file check.")
@click.pass_context
def format_kegg(ctx, kegg_pep_loc, gene_ko_link_loc, kegg_pep_root_dir, kegg_db, kegg_download_date, skip_gene_ko_link):
    """Format KEGG database for use in $workflow.manifest.name. Standalone operation, will exit after completion."""
    click.echo("Formatting KEGG database")

    if kegg_pep_loc and kegg_pep_root_dir:
        raise click.exceptions.BadArgumentUsage("Either include `kegg_pep_loc` or `kegg_pep_root_dir`, but not both.")
    if not kegg_pep_loc and not kegg_pep_root_dir:
        raise click.exceptions.BadArgumentUsage("Must include `kegg_pep_loc` or `kegg_pep_root_dir`.")

    output_dir = ctx.parent.params['output_dir']  # TODO, pass in outputdir for where to save kegg.pep
    threads = ctx.parent.params['threads']
    if kegg_pep_root_dir:

        p = subprocess.run(["bash", str(DRAM_DIR / "database/cat_kegg_pep.sh"), kegg_pep_root_dir, output_dir], capture_output=True)
        for line in p.stdout.decode().splitlines():
            click.echo(line)
        for line in p.stderr.decode().splitlines():
            click.echo(line, err=True)
        if p.returncode != 0:
            raise click.exceptions.ClickException("Error formatting KEGG database.")
        kegg_pep_loc = Path(output_dir) / "kegg.pep"
    format_kegg_database.prepare_databases(
        kegg_loc=kegg_pep_loc,
        gene_ko_link_loc=gene_ko_link_loc,
        output_dir=kegg_db,
        kegg_download_date=kegg_download_date,
        threads=threads,
        skip_gene_ko_link=skip_gene_ko_link
    )

if __name__ == '__main__':
    cli()