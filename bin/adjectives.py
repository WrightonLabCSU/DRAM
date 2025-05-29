#!/usr/bin/env python
"""Contains main entry point to the program and helper functions"""
import os
import click

import pandas as pd

from rule_adjectives.rule_graph import RuleParser, get_positive_genes
from rule_adjectives.annotations import Annotations


class PythonLiteralOption(click.Option):

    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)


def get_package_path(local_path):
    """
    Locate the package data or non python files

    :param local_path:
    :returns:
    """
    abs_snake_path = os.path.join(os.path.dirname(
        os.path.abspath(__file__)),
                                  "rule_adjectives",
        local_path)
    return abs_snake_path


def list_adjectives(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    rules = RuleParser(get_package_path('rules.tsv'), verbose=False)
    print("In the current rules file, these adjectives are available:")
    for i in rules.data.index[~rules.data['name'].isna()].unique():
        print(i)


def list_adjective_name(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    rules = RuleParser(get_package_path('rules.tsv'), verbose=False)
    print("In the current rules file, these adjectives are available:")
    for i in rules.data['name'].unique():
        print(i)

def show_rules_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(get_package_path('rules.tsv'))
    
def validate_comma_separated(ctx, param, value):
    print(value)
    if not value:
        return []
    return value.split(',')

@click.command()
@click.option('--annotations', type=click.Path(exists=True), required=True, help="One of only 2 required files. Path to a DRAM annotations file.")
@click.option('-o', '--output', type=click.Path(), default='adjectives.tsv', help="Path for the output table. A true false table created by this script.")
@click.option('-a', '--adjectives_list', default=[], callback=validate_comma_separated, help="A comma seperated list of adjectives ('adj1,adj2,adj3'), by name, to evaluate. This limits the number of adjectives that are evaluated and is faster.")
@click.option('-p', '--plot_adjectives', multiple=True, default=[], help="A list of adjectives, by name, to plot. This limits the number of adjectives that are plotted and is probably needed for speed.")
@click.option('-g', '--plot_genomes', multiple=True,
              default=[], )
@click.option('--plot_path', type=click.Path(exists=False),
              default=None,
              help='will become a folder of output plots, no path no plots.')
@click.option('--strainer_tsv', type=click.Path(exists=False), default=None, help='The path for a tsv that will pass to strainer to filter genes. The only option at this time is pgtb for positive genes that are on true bugs.')
@click.option('--strainer_type', type=click.Path(exists=False), default=None, help='The type of process that should make the strainer file.')
@click.option('--debug_ids_by_fasta_to_tsv', type=click.Path(exists=False), default=None,
              help='This is a tool to debug the list of IDs found by DRAM it is mostly for experts')
@click.option('--rules_tsv', type=click.Path(exists=True),
              default=get_package_path('rules.tsv'),
              help="This is an optional path to a rules file with strict formatting. It will over write the original rules file that is stored with the script.")
@click.option('--show_rules_path', is_flag=True, callback=show_rules_path,
              expose_value=False, is_eager=True,
              help="Show the path to the default rules path.")
@click.option('--list_name', is_flag=True, callback=list_adjective_name,
              expose_value=False, is_eager=True,
              help="List the names for all adjectives.tsv that are"
                   " available, you can pass these names to limit the"
                   " adjectives that are evaluated")
@click.option('--list_id', is_flag=True, callback=list_adjectives,
              expose_value=False, is_eager=True,
              help="List the names for all adjectives.tsv that are"
                   " available, you can pass these names to limit the"
                   " adjectives that are evaluated")
# @click.argument('-p', type=click.Path(exists=True))
def evaluate(annotations:str, output:str,
             rules_tsv:str=get_package_path('rules.tsv'),
             adjectives_list:list=None, plot_adjectives:list=None,
             plot_genomes:list=None,plot_path:str=None,
             debug_ids_by_fasta_to_tsv:str=None,
             strainer_tsv:str=None, strainer_type='pgtb'):
    """Using a DRAM annotations file make a table of adjectives."""
    rules = RuleParser(rules_tsv, verbose=False, adjectives=adjectives_list)
    annotations = Annotations(annotations)
    adjectives = rules.check_genomes(annotations)
    if debug_ids_by_fasta_to_tsv is not None:
        annotations.ids_by_fasta.to_csv(debug_ids_by_fasta_to_tsv, sep='\t')
        exit()
    adjectives.to_csv(output, sep='\t')
    # annotations.ids_by_fasta.iloc[1]['annotations']
    if plot_path is not None:
        rules.plot_cause(plot_path, adjectives=plot_adjectives,
                         genomes=plot_genomes, show_steps=False
                         )
    if strainer_tsv is not None:
        strainer_data = get_positive_genes(rules, annotations, adjectives)
        strainer_data.to_csv(strainer_tsv, sep='\t')


@click.command()
@click.argument('plot_path', type=click.Path(exists=False),
                default=None)#, help='will become a folder of output plots, no path no plots.')
@click.option('-a', '--adjectives', multiple=True, default=[], help="A list of adjectives, by name, to evaluate. This limits the number of adjectives that are evaluated and is faster.")
@click.option('--rules_tsv', type=click.Path(exists=True),
              default=get_package_path('rules.tsv'),
              help='The path that will become a folder of output plots, no path no plots.') # , help='The rules file which adhere to strict formating' )
@click.option('--list_name', is_flag=True, callback=list_adjective_name,
              expose_value=False, is_eager=True,
              help="List the names for all adjectives.tsv that are"
                   " available, you can pass these names to limit the"
                   " adjectives that are evaluated")
@click.option('--show_rules_path', is_flag=True, callback=show_rules_path,
              expose_value=False, is_eager=True,
              help="Show the path to the default rules path.")
@click.option('--list_id', is_flag=True, callback=list_adjectives,
              expose_value=False, is_eager=True,
              help="List the names for all adjectives.tsv that are"
                   " available, you can pass these names to limit the"
                   " adjectives that are evaluated")
# @click.argument('-p', type=click.Path(exists=True))
def rule_plot(rules_tsv:str=get_package_path('rules.tsv'),
              adjectives:list=None, plot_path:str=None):
    """
    Using a DRAM annotations file make a table of adjectives.
    """
    rules = RuleParser(rules_tsv, verbose=False, adjectives=adjectives)
    rules.plot_rule(plot_path)

if __name__ == "__main__":
    evaluate()