
import click
from functools import wraps

import sys

# # Import Python scripts as modules (example imports, replace with actual modules)
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


PIPELINE_ORDER = [
    'call', 'annotate', 'distill', 'product'
]

PIPELINES_WITH_FORMAT_KEGG = [
    'format_kegg', *PIPELINE_ORDER
]

# PIPELINE_ORDER = {
#     'format_kegg': [],
#     'call': [],
#     'annotate': [],
#     'distill': [],
#     'product': []
# }

@click.group(chain=True, invoke_without_command=True)
@click.pass_context
def cli(ctx):
    """Main CLI tool."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
    return

    ctx.ensure_object(dict)


def processor(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    wrapper.name = func.__name__
    return wrapper

@cli.command()
# @click.option('--call', is_flag=True, help="Call genes on the input FASTA files using Prodigal.")
@click.option('--input_fasta', type=click.Path(), default='./input_fasta/*.fa*',
              help="Directory containing input FASTA files. Default: './input_fasta/*.fa*'")
@click.option('--rename', is_flag=True, help=(
        "Rename FASTA headers based on file name. "
        "Example: sample1.fa --> (header renamed to) >sample1. "
        "If headers are already renamed, omit '--call'."
))
@click.option('--prodigal_mode', type=click.Choice(['single', 'meta']), default='single',
              help="Prodigal mode for gene calling. Default: 'single'.")
@click.option('--prodigal_tras_table', type=click.IntRange(1, 25), default=1,
              help="Translation table to use (1-25). Default: '1'.")
@click.option('--min_contig_len', type=int, default=2500,
              help="Minimum contig length in base pairs. Default: '2500'.")
@processor
def call(input_fasta, rename, prodigal_mode, prodigal_tras_table, min_contig_len):
    """Call genes using Prodigal."""
    click.echo("Calling genes")
    # TODO: Implement corresponding Python function


@cli.command()
# @click.option('--annotate', is_flag=True, help="Annotate called genes using downloaded databases.")
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
@click.option('--input_fasta', type=click.Path(), default='./input_fasta/',
              help="Directory containing input FASTA files. Default: './input_fasta/'.")
@click.option('--input_genes', type=click.Path(),
              help="Directory containing called genes (.faa).")
@click.option('--add_annotations', type=click.Path(),
              help="Path to old annotations to include in the current run.")
@click.option('--generate_gff', is_flag=True, help="Generate output GFF files for each sample.")
@click.option('--generate_gbk', is_flag=True, help="Generate output GBK files for each sample.")
@processor
def annotate(use_db, use_dbset, input_fasta, input_genes, add_annotations, generate_gff, generate_gbk):
    """Annotate called genes using downloaded databases"""
    click.echo("Annotating genes")
    # combine_annotations(input_dir=input_dir, output_dir=output_dir)


@cli.command()
# @click.option('--distill', is_flag=True, help="Distill annotations into a multi-sheet distillate.xlsx.")
@click.option('--annotations', type=click.Path(),
              help="Path to annotations.tsv file. Required if not running --call and --annotate.")
@click.option('--rrnas', type=click.Path(), help="Path to rRNA.tsv file.")
@click.option('--trnas', type=click.Path(), help="Path to tRNA.tsv file.")
@click.option('--bin_quality', type=click.Path(), help="Path to bin-quality.tsv file.")
@click.option('--taxa', type=click.Path(), help="Path to bin-taxonomy.tsv file.")
@click.option('--distill_topic', type=str, help=(
        "Topics for distillation (e.g., 'carbon transport'). "
        "If multiple topics, enclose in quotes."))
@click.option('--distill_ecosystem', type=str, help=(
        "Ecosystem for distillation (e.g., 'eng_sys ag'). "
        "If multiple ecosystems, enclose in quotes."))
@click.option('--distill_custom', type=click.Path(),
              help="Path to custom distillate TSV file.")
@processor
def distill(annotations, rrnas, trnas, bin_quality, taxa, distill_topic, distill_ecosystem, distill_custom):
    """Distill the annotations into a multi-sheet distillate.xlsx"""
    click.echo("Distilling annotations")
    # count_annotations(annotations=annotations, output=output)

@cli.command()
# @click.option('--product', is_flag=True, help="Generate a product visualization of annotations.")
@click.option('--annotations', type=click.Path(), required=True,
              help="Path to annotations.tsv file.")
@click.option('--groupby_column', type=str, default='fasta',
              help="Column for grouping in annotations file. Default: 'fasta'.")
@click.option('--module_steps_form', type=click.Path(),
              help="Path to module steps form TSV file.")
@click.option('--etc_module_form', type=click.Path(),
              help="Path to etc module form TSV file.")
@click.option('--function_heatmap_form', type=click.Path(),
              help="Path to function heatmap form TSV file.")
@processor
def product(annotations, groupby_column, module_steps_form, etc_module_form, function_heatmap_form):
    """Generate a product visualization of the annotations."""
    click.echo("Generating product visualization")

@cli.command()
# @click.option('--format_kegg', is_flag=True, help="Format KEGG database for use in the workflow.")
@click.option('--kegg_pep_loc', type=click.Path(), required=True,
              help="Path to KEGG peptide location.")
@click.option('--gene_ko_link_loc', type=click.Path(), help="Path to genes_ko_link file.")
@click.option('--kegg_db', type=click.Path(), default='./kegg_db',
              help="Path to save formatted KEGG database. Default: './kegg_db'.")
@click.option('--kegg_download_date', type=str, help="Date KEGG database was downloaded (YYYY-MM-DD).")
@click.option('--skip_gene_ko_link', is_flag=True, help="Skip the gene_ko_link file check.")
@processor
def format_kegg(kegg_pep_loc, gene_ko_link_loc, kegg_db, kegg_download_date, skip_gene_ko_link):
    """Format KEGG database for use in $workflow.manifest.name. Standalone operation, will exit after completion."""
    click.echo("Formatting KEGG database")

    # format_kegg_db(kegg_pep_loc=kegg_pep_loc, threads=threads)


if __name__ == '__main__':
    commands = sys.argv[1:]  # We can remove later because slicing creates a copy

    if "format-kegg" in commands and any(cmd in PIPELINE_ORDER for cmd in commands):
        raise click.exceptions.BadArgumentUsage("format_kegg must be run alone.")

    # we use sys.argv here instaed of commands since we are removing later in the loop,
    # and can't remove from the loop we are iterating from
    # We use set, to make sure we hit each cmd only once, so if it is cmd1 cmd2 cmd1, we only remove cmd1 once,
    # and detect the duplicate
    for cmd in set(sys.argv[1:]):
        if cmd in PIPELINE_ORDER:
            commands.remove(cmd)

    click.echo(commands)
    if commands:
        raise click.exceptions.BadArgumentUsage(f"Pipeline steps must be run in order: {', '.join(PIPELINES_WITH_FORMAT_KEGG)}")

    # If all is good, run it
    cli()
