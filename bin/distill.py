#!/usr/bin/env python
"""This is the script that distills the genomes"""
import logging
from collections import Counter
from itertools import chain
import pandas as pd
from collections import Counter, defaultdict
import re
import numpy as np
import click
import os
from pathlib import Path

from utils.logger import get_logger
from utils.click_utils import validate_comma_separated
from rule_adjectives.annotations import FUNCTION_DICT

# TODO: add RBH information to output
# TODO: add flag to output table and not xlsx
# TODO: add flag to output heatmap table

logger = get_logger()

COL_GENE_ID, COL_GENE_DESCRIPTION, COL_MODULE, COL_SHEET, COL_HEADER, COL_SUBHEADER = 'gene_id', 'gene_description', 'pathway', 'topic_ecosystem','category', 'subcategory'
FRAME_COLUMNS = [COL_GENE_ID, COL_GENE_DESCRIPTION, COL_MODULE, COL_SHEET, COL_HEADER, COL_SUBHEADER]
RRNA_TYPES = ['5S rRNA', '16S rRNA', '23S rRNA']
TAXONOMY_LEVELS = ['d', 'p', 'c', 'o', 'f', 'g', 's']
CONSTANT_DISTILLATE_COLUMNS = [COL_GENE_ID, COL_GENE_DESCRIPTION, COL_MODULE, COL_HEADER, COL_SUBHEADER]
DISTILATE_SORT_ORDER_COLUMNS = [COL_HEADER, COL_SUBHEADER, COL_MODULE, COL_GENE_ID]
EXCEL_MAX_CELL_SIZE = 32767

FASTA_COLUMN = os.getenv('FASTA_COLUMN')
DISTILL_DIR = Path(__file__).parent.parent / "assets/forms/distill_sheets"


def check_columns(data, logger):
    functions = {i:j for i,j in FUNCTION_DICT.items() if i in data.columns}
    missing = [i for i in FUNCTION_DICT if i not in data.columns]
    logger.info("Note: the following id fields "
          f"were not in the annotations file and are not being used: {missing},"
          f" but these are {list(functions.keys())}")


def get_ids_from_annotations_by_row(data):
    functions = {i:j for i,j in FUNCTION_DICT.items() if i in data.columns}
    out = data.apply(lambda x: {i for k, v in functions.items() if not pd.isna(x[k])
                          for i in v(str(x[k])) if not pd.isna(i)}, axis=1)
    return out


def get_ids_from_annotations_all(data):
    data =  get_ids_from_annotations_by_row(data)
    data.apply(list)
    out = Counter(chain(*data.values))
    return out


def fill_genome_summary_frame(annotations, genome_summary_frame, groupby_column, logger):
    genome_summary_id_sets = [set([str(k).strip() for k in j.split(',')]) for j in genome_summary_frame[COL_GENE_ID]]
    logger.info(f"Genome summary ID sets: {genome_summary_id_sets}")
    
    def fill_a_frame(frame: pd.DataFrame):
        id_dict = get_ids_from_annotations_all(frame)
        logger.info(f"ID dictionary for {frame.name}: {id_dict}")
        
        counts = list()
        for set_ in genome_summary_id_sets:
            identifier_count = 0
            for gene_id in set_:
                # Try matching with and without '.hmm'
                matching_keys = [key for key in id_dict.keys() if gene_id == key or key.startswith(gene_id + ".")]
                for key in matching_keys:
                    identifier_count += id_dict[key]
            counts.append(identifier_count)
        # logger.info(f"Counts for {frame.name}: {counts}")
        
        return pd.Series(counts, index=genome_summary_frame.index)
    
    counts = annotations.groupby(groupby_column, sort=False).apply(fill_a_frame)
    genome_summary_frame = pd.concat([genome_summary_frame, counts.T], axis=1)
    
    return genome_summary_frame


def fill_genome_summary_frame_gene_names(annotations, genome_summary_frame, groupby_column, logger):
    genome_summary_id_sets = [set([k.strip() for k in j.split(',')]) for j in genome_summary_frame[COL_GENE_ID]]
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        # make dict of identifiers to gene names
        id_gene_dict = defaultdict(list)
        for gene, ids in get_ids_from_annotations_by_row(frame).items():
            for id_ in ids:
                id_gene_dict[id_].append(gene)
        # fill in genome summary_frame
        values = list()
        for id_set in genome_summary_id_sets:
            this_value = list()
            for id_ in id_set:
                this_value += id_gene_dict[id_]
            values.append(','.join(this_value))
        genome_summary_frame[genome] = values
    return genome_summary_frame


def summarize_rrnas(rrnas_df, groupby_column=FASTA_COLUMN):
    genome_rrna_dict = dict()
    for genome, frame in rrnas_df.groupby(groupby_column):
        genome_rrna_dict[genome] = Counter(frame['type'])
    row_list = list()
    for rna_type in RRNA_TYPES:
        row = [rna_type, '%s ribosomal RNA gene' % rna_type.split()[0], 'rRNA', 'rRNA', '', '']
        for genome, rrna_dict in genome_rrna_dict.items():
            row.append(genome_rrna_dict[genome].get(rna_type, 0))
        row_list.append(row)
    rrna_frame = pd.DataFrame(row_list, columns=FRAME_COLUMNS + list(genome_rrna_dict.keys()))
    return rrna_frame


def summarize_trnas(trnas_df, groupby_column=FASTA_COLUMN):
    # first build the frame
    combos = {(line.type, line.codon, line.note) for _, line in trnas_df.iterrows()}
    frame_rows = list()
    for combo in combos:
        if combo[2] == 'pseudo':
            gene_id = '%s, pseudo (%s)'
            gene_description = '%s pseudo tRNA with %s Codon'
        else:
            gene_id = '%s (%s)'
            gene_description = '%s tRNA with %s Codon'
        gene_id = gene_id % (combo[0], combo[1])
        gene_description = gene_description % (combo[0], combo[1])
        module_description = '%s tRNA' % combo[0]
        frame_rows.append([gene_id, gene_description, module_description, 'tRNA', 'tRNA', ''])
    trna_frame = pd.DataFrame(frame_rows, columns=FRAME_COLUMNS)
    trna_frame = trna_frame.sort_values(COL_GENE_ID)
    # then fill it in
    trna_frame = trna_frame.set_index(COL_GENE_ID)
    for group, frame in trnas_df.groupby(groupby_column):
        gene_ids = list()
        for index, line in frame.iterrows():
            if line.note == 'pseudo':
                gene_id = '%s, pseudo (%s)'
            else:
                gene_id = '%s (%s)'
            gene_ids.append(gene_id % (line.type, line.codon))
        trna_frame[group] = pd.Series(Counter(gene_ids))
    trna_frame = trna_frame.reset_index()
    trna_frame = trna_frame.fillna(0)
    return trna_frame


def make_genome_summary(annotations, genome_summary_frame, logger, trna_frame=None, rrna_frame=None, groupby_column=FASTA_COLUMN):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(fill_genome_summary_frame(annotations, genome_summary_frame.copy(), groupby_column, logger))

    # add rRNAs
    if rrna_frame is not None:
        summary_frames.append(summarize_rrnas(rrna_frame, groupby_column))

    # add tRNAs
    if trna_frame is not None:
        summary_frames.append(summarize_trnas(trna_frame, groupby_column))

    # merge summary frames
    summarized_genomes = pd.concat(summary_frames, sort=False)
    return summarized_genomes


def split_column_str(names):
    if len(names) < EXCEL_MAX_CELL_SIZE:
        return [names]
    out = ['']
    name_list = names.split(',')
    j = 0
    for i in name_list:
        if len(out[j]) + len(i) + 1 < EXCEL_MAX_CELL_SIZE:
            out[j] = ','.join([out[j], i])
        else:
            j += 1
            out += ['']
    return out


def split_names_to_long(col:pd.Series):
    dex = col.index
    splits = [split_column_str(i) for i in col.values]
    ncols =  max([len(i) for i in splits])
    col_names = [col.name if i == 0 else f"{col.name}[{i + 1}]" for i in range(ncols)]
    return pd.DataFrame(splits, columns=col_names, index=dex).fillna('')


def write_summarized_genomes_to_xlsx(summarized_genomes, output_file):
    # turn all this into an xlsx
    with pd.ExcelWriter(output_file) as writer:
        for sheet, frame in summarized_genomes.groupby(COL_SHEET, sort=False):
            frame = frame.sort_values(DISTILATE_SORT_ORDER_COLUMNS)
            frame = frame.drop([COL_SHEET], axis=1)
            gene_columns = list(set(frame.columns) - set(CONSTANT_DISTILLATE_COLUMNS))
            split_genes = pd.concat([split_names_to_long(frame[i].astype(str)) for i in gene_columns], axis=1)
            frame = pd.concat([frame[CONSTANT_DISTILLATE_COLUMNS],  split_genes], axis=1)
            frame.to_excel(writer, sheet_name=sheet, index=False)


@click.command()
@click.option("-i", "--input_file", required=True, help="Annotations path")
# @click.option("-o", "--output_dir", required=True, help="Directory to write summarized genomes")
@click.option('--log_file_path', 
                                    help="A name and loctation for the log file")
@click.option("--rrna_path", help="rRNA output from annotation")
@click.option("--trna_path", help="tRNA output from annotation")
@click.option("--groupby_column", help="Column from annotations to group as organism units",
                            default=FASTA_COLUMN)
@click.option("--distil_topics", default="default", help="Default distillates topics to run.")
@click.option("--distil_ecosystem", default="eng_sys,ag", help="Default distillates ecosystems to run.")
@click.option("--custom_distillate", default=[], callback=validate_comma_separated, help="Custom distillate forms to add your own modules, comma separated. ")
@click.option("--distillate_gene_names", is_flag=True,
    show_default=True, default=False,
                            help="Give names of genes instead of counts in genome metabolism summary")
def distill(input_file, trna_path=None, rrna_path=None, distil_topics=None, distil_ecosystem=None,
                      groupby_column=FASTA_COLUMN, log_file_path=None, custom_distillate=None,
                      distillate_gene_names=False):
    """Summarize metabolic content of annotated genomes"""
    # make output folder
    # mkdir(output_dir)

    # read in data
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    if 'bin_taxnomy' in annotations:
        annotations = annotations.sort_values('bin_taxonomy')

    # Check the columns are present
    check_columns(annotations, logger)

    if trna_path is None:
        trna_frame = None
    else:
        trna_frame = pd.read_csv(trna_path, sep='\t')
    if rrna_path is None:
        rrna_frame = None
    else:
        rrna_frame = pd.read_csv(rrna_path, sep='\t')


    distil_dir = Path(__file__).parent.parent / "assets/forms/distill_sheets"
    distil_sheets_names = []
    if "default" in distil_topics:
        distil_sheets_names = [
            DISTILL_DIR / "distill_carbon.tsv",
            DISTILL_DIR / "distill_energy.tsv",
            DISTILL_DIR / "distill_misc.tsv",
            DISTILL_DIR / "distill_nitrogen.tsv",
            DISTILL_DIR / "distill_transport.tsv"
        ]
    else:
        if 'carbon' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_carbon.tsv")
        if 'energy' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_energy.tsv")
        if 'misc' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_misc.tsv")
        if 'nitrogen' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_nitrogen.tsv")
        if 'transport' in distil_topics:
            distil_sheets_names.append(DISTILL_DIR / "distill_transport.tsv")
    
        
    if "ag" in distil_ecosystem:
        distil_sheets_names.append(DISTILL_DIR / "distill_ag.tsv")
    if "eng_sys" in distil_ecosystem:
        distil_sheets_names.append(DISTILL_DIR / "distill_eng_sys.tsv")
    
    if "camper_id" in annotations and ("default" in distil_topics or "camper" in distil_topics):
        distil_sheets_names.append(DISTILL_DIR / "distill_camper.tsv")
        
        
    if custom_distillate:
        for custom_sheet in custom_distillate:
            distil_sheets_names.append(custom_sheet)
    
    genome_summary_form = pd.concat(
        [pd.read_csv(sheet, 
                     sep='\t', 
                     usecols=FRAME_COLUMNS) 
         for sheet in distil_sheets_names],
    )
    
    logger.info('Retrieved distillate genome summary form')

    genome_summary_form = genome_summary_form.reset_index(drop=True)

    # make genome metabolism summary
    genome_summary = 'distillate.xlsx'
    if distillate_gene_names:
        logger.info(f'distillate_gene_names flag is {distillate_gene_names}. Giving gene names instead of counts in genome metabolism summary')
        summarized_genomes = fill_genome_summary_frame_gene_names(annotations, genome_summary_form, groupby_column, logger)
    else:
        logger.info(f'distillate_gene_names flag is {distillate_gene_names}. Giving counts instead of gene names in genome metabolism summary')
        summarized_genomes = make_genome_summary(annotations, genome_summary_form, logger, trna_frame, rrna_frame,
                                                 groupby_column)
    summarized_genomes.to_csv('summarized_genomes.tsv', sep='\t', index=None)
    write_summarized_genomes_to_xlsx(summarized_genomes, genome_summary)
    logger.info('Generated genome metabolism summary')

    
if __name__ == "__main__":
    distill()