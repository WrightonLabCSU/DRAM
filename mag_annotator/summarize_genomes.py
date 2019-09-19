import pandas as pd
from collections import Counter, defaultdict
from os import path, mkdir
import re
import altair as alt
import networkx as nx

from mag_annotator.utils import get_database_locs

# TODO: add RBH information to output
# TODO: add measure of redendancy of genes
# TODO: add total number of copies
# TODO: add tqdm progress bar
# TODO: add ability to take in GTDBTK file and add taxonomy to annotations

FRAME_COLUMNS = ['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader']


def get_all_ids(frame):
    id_counts = dict()
    # kegg KO ids
    id_counts.update(Counter([j for i in frame.kegg_id if not pd.isna(i) for j in i.split(',')]))
    # cazy ids
    id_counts.update(Counter([j.strip().split('_')[0]
                              for i in frame.cazy_hits if not pd.isna(i) for j in i.split(';')]))
    # peptidase ids
    id_counts.update(Counter([i for i in frame.peptidase_family if not pd.isna(i)]))
    return id_counts


def fill_genome_summary_frame(annotations, group_column, genome_summary_frame):
    for genome, frame in annotations.groupby(group_column, sort=False):
        id_dict = get_all_ids(frame)
        genome_summary_frame[genome] = [id_dict[i] if i in id_dict else 0 for i in genome_summary_frame.gene_id]
    return genome_summary_frame


RRNA_TYPES = ['5S rRNA', '16S rRNA', '23S rRNA']


def summarize_rrnas(rrnas_df, groupby_column='fasta'):
    genome_rrna_dict = dict()
    for genome, frame in rrnas_df.groupby(groupby_column):
        genome_rrna_dict[genome] = Counter(frame['type'])
    row_list = list()
    for rna_type in RRNA_TYPES:
        row = [rna_type, '%s ribosomal RNA gene' % rna_type.split()[0], 'rRNA', 'ribosomal RNA genes', 'rRNA genes',
               True]
        for genome, rrna_dict in genome_rrna_dict.items():
            row.append(genome_rrna_dict[genome].get(type, 0))
        row_list.append(row)
    rrna_frame = pd.DataFrame(row_list, columns=FRAME_COLUMNS+list(genome_rrna_dict.keys()))
    return rrna_frame


def summarize_trnas(trnas_df, groupby_column='fasta'):
    # first build the frame
    combos = set()
    for index, line in trnas_df.iterrows():
        combos.add((line.Type, line.Codon, line.Note))
    frame_rows = list()
    for combo in combos:
        if combo[2] == 'pseudo':
            gene_id = '%s, pseudo (%s)'
            gene_description = '%s pseudo tRNA with %s Codon'
        else:
            gene_id = '%s (%s)'
            gene_description = '%s pseudo tRNA with %s Codon'
        gene_id = gene_id % (combo[0], combo[1])
        gene_description = gene_description % (combo[0], combo[1])
        module_id = combo[0]
        module_description = '%s tRNA' % combo[0]
        module_family = 'tRNA genes'
        frame_rows.append([gene_id, gene_description, module_id, module_description, module_family, ''])
    trna_frame = pd.DataFrame(frame_rows, columns=FRAME_COLUMNS)
    trna_frame = trna_frame.sort_values('gene_id')
    # then fill it in
    trna_frame = trna_frame.set_index('gene_id')
    for group, frame in trnas_df.groupby(groupby_column):
        gene_ids = list()
        for index, line in frame.iterrows():
            if line.Note == 'pseudo':
                gene_id = '%s, pseudo (%s)'
            else:
                gene_id = '%s (%s)'
            gene_ids.append(gene_id % (line.Type, line.Codon))
        trna_frame[group] = pd.Series(Counter(gene_ids))
    trna_frame = trna_frame.reset_index()
    trna_frame = trna_frame.fillna(0)
    return trna_frame


def make_genome_summary(annotations, genome_summary_frame, trna_frame=None, rrna_frame=None,
                        group_column='fasta', viral=False):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(fill_genome_summary_frame(annotations, group_column, genome_summary_frame.copy()))

    # add rRNAs
    if rrna_frame is not None:
        summary_frames.append(summarize_rrnas(rrna_frame, group_column))

    # add tRNAs
    if trna_frame is not None:
        summary_frames.append(summarize_trnas(trna_frame, group_column))

    # merge summary frames
    summarized_genomes = pd.concat(summary_frames, sort=False)

    # post processing
    if viral:  # filter out empty rows and columns if viral
        summarized_genomes_numbers_only = summarized_genomes[summarized_genomes.columns[7:]]
        # remove all zero rows for viral
        summarized_genomes = summarized_genomes.loc[summarized_genomes_numbers_only.sum(axis=1) > 0]
        # remove all zero columns so viruses with no AMGs
        good_columns = summarized_genomes_numbers_only.columns[summarized_genomes_numbers_only.sum(axis=0) > 0]
        summarized_genomes = summarized_genomes[list(summarized_genomes.columns[:7]) + list(good_columns)]

    return summarized_genomes


def make_genome_stats(annotations, rrna_frame=None, trna_frame=None, group_column='fasta'):
    rows = list()
    columns = ['genome', 'number of scaffolds', 'taxonomy', 'completeness', 'contamination']
    if rrna_frame is not None:
        columns += RRNA_TYPES
    if trna_frame is not None:
        columns.append('tRNA count')
    for genome, frame in annotations.groupby(group_column, sort=False):
        row = [genome, len(set(frame['scaffold']))]
        if 'bin_taxonomy' in frame:
            row.append(frame['bin_taxonomy'][0])
        if 'bin_completeness' in frame:
            row.append(frame['bin_completeness'][0])
        if 'bin_contamination' in frame:
            row.append(frame['bin_contamination'][0])
        if rrna_frame is not None:
            genome_rrnas = rrna_frame.loc[rrna_frame.fasta == genome]
            for rrna in RRNA_TYPES:
                sixteens = genome_rrnas.loc[genome_rrnas.type == rrna]
                if sixteens.shape[0] == 0:
                    row.append('')
                elif sixteens.shape[0] == 1:
                    row.append('%s, (%s, %s)' % (sixteens.index[0], sixteens.begin[0], sixteens.end[0]))
                else:
                    row.append('%s present' % sixteens.shape[0])
        if trna_frame is not None:
            row.append(trna_frame.loc[trna_frame[group_column] == genome].shape[0])
        rows.append(row)
    genome_stats = pd.DataFrame(rows, columns=columns)
    return genome_stats


def get_ids_from_annotation(frame):
    id_list = list()
    # get kegg ids
    id_list += [j for i in frame.kegg_id.dropna() for j in i.split(',')]
    # get ec numbers
    for kegg_hit in frame.kegg_hit.dropna():
        id_list += [i[1:-1] for i in re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', kegg_hit)]
    # get merops ids
    id_list += [j for i in frame.peptidase_family.dropna() for j in i.split(';')]
    # get cazy ids
    id_list += [j.split(' ')[0] for i in frame.cazy_hits.dropna() for j in i.split(';')]
    return id_list


def build_module_net(module_df):
    # build net from a set of module paths
    num_steps = max([int(i.split(',')[0]) for i in set(module_df.path)]) + 1
    module_net = nx.DiGraph(num_steps=num_steps, module_id=list(module_df.module)[0],
                            module_name=list(module_df.module_name)[0])
    for module_path, frame in module_df.groupby('path'):
        split_path = [int(i) for i in module_path.split(',')]
        module_net.add_node(module_path, kos=set(frame.ko))
        # add incoming edge
        if module_path[0] == 0:
            module_net.add_edge('begin', module_path)
        else:
            module_net.add_edge('end_step_%s' % (split_path[0]-1), module_path)
        # add outgoing edge
        if split_path[0] == num_steps:
            module_net.add_edge(module_path, 'end')
        else:
            module_net.add_edge(module_path, 'end_step_%s' % split_path[0])
    return module_net


def get_module_coverage(kos, module_net):
    # prune network based on what kos were observed
    pruned_module_net = module_net.copy()
    module_kos_present = set()
    for node, data in module_net.nodes.items():
        if 'kos' in data:
            ko_overlap = data['kos'] & kos
            if len(ko_overlap) == 0:
                pruned_module_net.remove_node(node)
            else:
                module_kos_present = module_kos_present | ko_overlap
    # count number of missing steps
    missing_steps = list()
    for node, data in pruned_module_net.nodes.items():
        if ('end_step' in node) and pruned_module_net.in_degree(node) == 0:
            missing_steps.append(int(node.split('_')[-1]))
    # get statistics
    num_steps = pruned_module_net.graph['num_steps']
    num_steps_present = num_steps-len(missing_steps)
    coverage = num_steps_present/num_steps
    return num_steps, num_steps_present, coverage, sorted(module_kos_present)


def make_module_coverage_df(annotation_df, module_nets):
    kos_to_genes = defaultdict(list)
    for gene_id, ko_list in annotation_df.kegg_id.iteritems():
        if type(ko_list) is str:
            for ko in ko_list.split(','):
                kos_to_genes[ko].append(gene_id)
    coverage_dict = {}
    for i, (module, net) in enumerate(module_nets.items()):
        module_steps, module_steps_present, module_coverage, module_kos = get_module_coverage(set(kos_to_genes.keys()),
                                                                                              net)
        module_genes = sorted([gene for ko in module_kos for gene in kos_to_genes[ko]])
        coverage_dict[module] = [net.graph['module_name'], module_steps, module_steps_present, module_coverage,
                                 len(module_kos), ','.join(module_kos), ','.join(module_genes)]
    coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index',
                                         columns=['module_name', 'steps', 'steps_present', 'step_coverage', 'ko_count',
                                                  'kos_present', 'genes_present'])
    return coverage_df


HEATMAP_CELL_HEIGHT = 10
HEATMAP_CELL_WIDTH = 10


def make_module_coverage_heatmap(annotations, module_nets, mag_order=None, groupby_column='fasta'):
    # go through each scaffold to check for modules
    module_coverage_dict = dict()
    for scaffold, frame in annotations.groupby(groupby_column, sort=False):
        module_coverage_dict[scaffold] = make_module_coverage_df(frame, module_nets)
    module_coverage = pd.concat(module_coverage_dict)
    num_mags_in_frame = len(set(annotations[groupby_column]))

    c = alt.Chart(module_coverage).encode(
        x=alt.X('MAG'),
        y=alt.Y('module_name', title='Module', sort=mag_order),
        tooltip=[alt.Tooltip('MAG', title='MAG'),
                 alt.Tooltip('module_name', title='Module Name'),
                 alt.Tooltip('steps', title='Module steps'),
                 alt.Tooltip('steps_present', title='Steps present')
                 ]
    )

    module_coverage_heatmap = c.mark_rect().encode(color='step_coverage').properties(
        width=HEATMAP_CELL_WIDTH * num_mags_in_frame,
        height=HEATMAP_CELL_HEIGHT * len(module_nets))

    return module_coverage_heatmap


def get_ordered_uniques(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def make_functional_heatmap(annotations, function_heatmap_form, groupby_column='fasta'):
    # clean up function heatmap form
    function_heatmap_form = function_heatmap_form.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    function_heatmap_form = function_heatmap_form.fillna('')
    # build dict of ids per genome
    genome_to_id_dict = dict()
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        id_list = get_ids_from_annotation(frame)
        genome_to_id_dict[genome] = set(id_list)
    # build long from data frame
    rows = list()
    for _, row in function_heatmap_form.iterrows():
        function_id_set = set([i.strip() for i in row.function_ids.strip().split(',')])
        for bin_name, id_set in genome_to_id_dict.items():
            functions_present = set.intersection(id_set, function_id_set)
            present_in_bin = len(functions_present) > 0
            rows.append([row.category, row.subcategory, row.function_name, ', '.join(functions_present),
                         row.long_function_name, row.gene_symbol, bin_name, present_in_bin])
    long_frame = pd.DataFrame(rows, columns=list(function_heatmap_form.columns) + ['bin', 'present'])
    # build heatmaps
    charts = list()
    grouped_function_names = long_frame.groupby('category', sort=False)
    mag_order = get_ordered_uniques(annotations.sort_values('bin_taxonomy')['fasta'])
    for i, (group, frame) in enumerate(grouped_function_names):
        # set variables for chart
        num_function_names_in_category = len(set(frame.function_name))
        num_mags_in_frame = len(set(frame.bin))
        chart_width = HEATMAP_CELL_WIDTH * num_function_names_in_category
        chart_height = HEATMAP_CELL_HEIGHT * num_mags_in_frame
        function_order = list(function_heatmap_form.loc[function_heatmap_form['category'] == group].function_name)
        # if this is the first chart then make y-ticks otherwise none
        if i == 0:
            y = alt.Y('bin', title=None, axis=alt.Axis(labelLimit=0), sort=mag_order)
        else:
            y=alt.Y('bin', axis=alt.Axis(title=None, labels=False, ticks=False), sort=mag_order)
        # set up colors for chart
        rect_colors = alt.Color('present',
                                legend=alt.Legend(title="Function is Present",
                                symbolType='square',
                                values=[True, False]),
                                sort=[True, False],
                                scale=alt.Scale(range=['#e5f5f9', '#2ca25f']))
        # define chart
        c = alt.Chart(frame, title=alt.TitleParams(group)).encode(  # TODO: Figure out how to angle title
            x=alt.X('function_name', title=None, axis=alt.Axis(labelLimit=0, labelAngle=90), sort=function_order),
            tooltip=[alt.Tooltip('bin', title='MAG'),
                     alt.Tooltip('category', title='Category'),
                     alt.Tooltip('subcategory', title='Subcategory'),
                     alt.Tooltip('function_ids', title='Function IDs'),
                     alt.Tooltip('function_name', title='Function'),
                     alt.Tooltip('long_function_name', title='Description'),
                     alt.Tooltip('gene_symbol', title='Gene Symbol')]
        ).mark_rect().encode(y=y, color=rect_colors).properties(
            width=chart_width,
            height=chart_height)
        charts.append(c)
    # merge and return
    function_heatmap = alt.hconcat(*charts)
    return function_heatmap


HEATMAP_MODULES = ['M00001', 'M00009', 'M00004']


def summarize_genomes(input_file, trna_path, rrna_path, output_dir, groupby_column, viral=False):
    # read in data
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    if 'bin_taxnomy' in annotations:
        annotations = annotations.sort_values('bin_taxonomy')

    if trna_path is None:
        trna_frame = None
    else:
        trna_frame = pd.read_csv(trna_path, sep='\t', index_col=0)
    if rrna_path is None:
        rrna_frame = None
    else:
        rrna_frame = pd.read_csv(rrna_path, sep='\t', index_col=0)

    # get db_locs and read in dbs
    db_locs = get_database_locs()
    if 'genome_summary_form' not in db_locs:
        raise ValueError('Genome summary form location must be set in order to summarize genomes')
    if 'module_step_form' not in db_locs:
        raise ValueError('Module step form location must be set in order to summarize genomes')
    if 'function_heatmap_form' not in db_locs:
        raise ValueError('Functional heat map form location must be set in order to summarize genomes')

    # read in dbs
    genome_summary_form = pd.read_csv(db_locs['genome_summary_form'], sep='\t')
    module_steps_form = pd.read_csv(db_locs['module_step_form'], sep='\t')
    function_heatmap_form = pd.read_csv(db_locs['function_heatmap_form'], sep='\t')

    # make output folder
    mkdir(output_dir)

    # make genome stats
    if not viral:
        genome_stats = make_genome_stats(annotations, rrna_frame, trna_frame, groupby_column)
        genome_stats.to_csv(path.join(output_dir, 'genome_stats.tsv'), sep='\t', index=False)

    # make genome metabolism summary
    genome_summary = make_genome_summary(annotations, genome_summary_form, trna_frame, rrna_frame, groupby_column,
                                         viral)
    genome_summary.to_csv(path.join(output_dir, 'genome_metabolism_summary.tsv'), sep='\t', index=False)

    # make heatmaps
    module_nets = {module: build_module_net(module_df)
                   for module, module_df in module_steps_form.groupby('module') if module in HEATMAP_MODULES}
    module_coverage_heatmap = make_module_coverage_heatmap(annotations, module_nets, groupby_column)

    function_heatmap = make_functional_heatmap(annotations, function_heatmap_form, groupby_column)

    alt.vconcat(module_coverage_heatmap, function_heatmap).save('heatmap.html')
