import pandas as pd
import altair as alt
import re
from os import path, mkdir
from functools import partial
from collections import defaultdict

from mag_annotator.utils import get_database_locs
from mag_annotator.summarize_genomes import get_ordered_uniques

VIRAL_DISTILLATE_COLUMNS = ['gene', 'scaffold', 'gene_id', 'gene_description', 'category', 'header',
                            'subheader', 'module', 'auxiliary_score', 'amg_flags']
VIRAL_LIQUOR_HEADERS = ['Category', 'Function', 'AMG Genes', 'Genes Present', 'VGF Name', 'Present in VGF']
HEATMAP_CELL_HEIGHT = 10
HEATMAP_CELL_WIDTH = 10

defaultdict_list = partial(defaultdict, list)


def filter_to_amgs(annotations, max_aux=4, remove_transposons=True, remove_fs=False):
    potential_amgs = list()
    for gene, row in annotations.iterrows():
        amg_flags = row['amg_flags']
        if not pd.isna(amg_flags):
            if ('V' not in amg_flags) and ('M' in amg_flags) and \
               (row['auxiliary_score'] <= max_aux) and ('A' not in amg_flags):
                if (remove_transposons and 'T' not in amg_flags) or not remove_transposons:
                    if (remove_fs and 'F' not in amg_flags) or not remove_fs:
                        potential_amgs.append(gene)
    return annotations.loc[potential_amgs]


def get_ids_from_row(row):
    id_list = list()
    # get kegg ids
    if not pd.isna(row.kegg_id):
        id_list += [j for j in row.kegg_id.split(',')]
    # get ec numbers
    if not pd.isna(row.kegg_hit):
        id_list += [i[1:-1] for i in re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', row.kegg_hit)]
    # get merops ids
    if not pd.isna(row.peptidase_family):
        id_list += [j for j in row.peptidase_family.split(';')]
    # get cazy ids
    if not pd.isna(row.cazy_hits):
        id_list += [j.strip().split(' ')[0] for j in row.cazy_hits.split(';')]
    return set(id_list)


def make_viral_distillate(potential_amgs, genome_summary_frame):
    rows = list()
    for gene, row in potential_amgs.iterrows():
        gene_ids = get_ids_from_row(row) & set(genome_summary_frame.index)
        for gene_id in gene_ids:
            gene_summary = genome_summary_frame.loc[gene_id]
            if type(gene_summary) is pd.Series:
                rows.append([gene, row['scaffold'], gene_id, gene_summary['gene_description'], gene_summary['sheet'],
                             gene_summary['header'], gene_summary['subheader'], gene_summary['module'],
                             row['auxiliary_score'], row['amg_flags']])
            else:
                for sub_gene_id, sub_gene_summary in genome_summary_frame.loc[gene_id].iterrows():
                    rows.append([gene, row['scaffold'], gene_id, sub_gene_summary['gene_description'],
                                 sub_gene_summary['sheet'], sub_gene_summary['header'],
                                 sub_gene_summary['subheader'], sub_gene_summary['module'],
                                 row['auxiliary_score'], row['amg_flags']])
    return pd.DataFrame(rows, columns=VIRAL_DISTILLATE_COLUMNS)


def make_viral_functional_df(annotations, genome_summary_frame, groupby_column='scaffold'):
    # build dict of ids per genome
    vgf_to_id_dict = defaultdict(defaultdict_list)
    for vgf, frame in annotations.groupby(groupby_column, sort=False):
        for gene, row in frame.iterrows():
            id_list = get_ids_from_row(row)
            for id_ in id_list:
                vgf_to_id_dict[vgf][id_].append(gene)
    # build long from data frame
    rows = list()
    for category, category_frame in genome_summary_frame.groupby('sheet'):
        for header, header_frame in category_frame.groupby('module'):
            header_id_set = set(header_frame.index.to_list())
            curr_rows = list()
            for vgf, id_dict in vgf_to_id_dict.items():
                present_in_bin = False
                functions_present = list()
                amgs_present = list()
                for id_, amgs in id_dict.items():
                    if id_ in header_id_set:
                        present_in_bin = True
                        functions_present.append(id_)
                        amgs_present += amgs
                curr_rows.append([category, header, ', '.join(amgs_present), ', '.join(functions_present), vgf, present_in_bin])
            if sum([i[-1] for i in curr_rows]) > 0:
                rows += curr_rows
    return pd.DataFrame(rows, columns=VIRAL_LIQUOR_HEADERS)


def make_viral_functional_heatmap(functional_df):
    # build heatmaps
    charts = list()
    for i, (group, frame) in enumerate(functional_df.groupby('Category', sort=False)):
        # set variables for chart
        function_order = get_ordered_uniques(list(frame['Function']))
        num_vgfs_in_frame = len(set(frame['VGF Name']))
        chart_width = HEATMAP_CELL_WIDTH * len(function_order)
        chart_height = HEATMAP_CELL_HEIGHT * num_vgfs_in_frame
        # if this is the first chart then make y-ticks otherwise none
        if i == 0:
            y = alt.Y('VGF Name', title=None, axis=alt.Axis(labelLimit=0))
        else:
            y = alt.Y('VGF Name', axis=alt.Axis(title=None, labels=False, ticks=False))
        # set up colors for chart
        rect_colors = alt.Color('Present in VGF',
                                legend=alt.Legend(symbolType='square', values=[True, False]),
                                sort=[True, False],
                                scale=alt.Scale(range=['#e5f5f9', '#2ca25f']))
        # define chart
        # TODO: Figure out how to angle title to take up less space
        c = alt.Chart(frame, title=alt.TitleParams(group)).encode(
            x=alt.X('Function', title=None, axis=alt.Axis(labelLimit=0, labelAngle=90), sort=function_order),
            tooltip=[alt.Tooltip('VGF Name'),
                     alt.Tooltip('Category'),
                     alt.Tooltip('Function'),
                     alt.Tooltip('AMG Genes'),
                     alt.Tooltip('Genes Present')]
        ).mark_rect().encode(y=y, color=rect_colors).properties(
            width=chart_width,
            height=chart_height)
        charts.append(c)
    # merge and return
    function_heatmap = alt.hconcat(*charts)
    return function_heatmap


def summarize_vgfs(input_file, output_dir, groupby_column='scaffold', max_auxiliary_score=4, remove_transposons=True,
                   remove_fs=False):
    # set up
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    annotations = annotations.fillna('')
    db_locs = get_database_locs()
    if 'genome_summary_form' not in db_locs:
        raise ValueError('Genome summary form location must be set in order to summarize genomes')
    mkdir(output_dir)
    genome_summary_form = pd.read_csv(db_locs['genome_summary_form'], sep='\t', index_col=0)
    # get potential AMGs
    potential_amgs = filter_to_amgs(annotations, max_aux=max_auxiliary_score, remove_transposons=remove_transposons,
                                    remove_fs=remove_fs)
    # make distillate
    viral_distillate = make_viral_distillate(potential_amgs, genome_summary_form)
    viral_distillate.to_csv(path.join(output_dir, 'vgf_amg_summary.tsv'), sep='\t', index=None)
    # make liquor
    viral_function_df = make_viral_functional_df(potential_amgs, genome_summary_form,
                                                 groupby_column=groupby_column)
    viral_functional_heatmap = make_viral_functional_heatmap(viral_function_df)
    viral_functional_heatmap.save(path.join(output_dir, 'vgf_amg_heatmap.html'))
