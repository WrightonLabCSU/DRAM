import pandas as pd
from collections import Counter, defaultdict
import networkx as nx
from os import path, mkdir

from mag_annotator.utils import get_database_locs

# TODO: add RBH information to output
# TODO: add measure of redendancy of genes
# TODO: add total number of copies
# TODO: add tqdm progress bar
# TODO: add ability to take in GTDBTK file and add taxonomy to annotations

FRAME_COLUMNS = ['gene_id', 'gene_description', 'module_id', 'module_description', 'module_family',
                 'key_gene']


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
    grouped = annotations.groupby(group_column)
    for genome, frame in grouped:
        id_dict = get_all_ids(frame)
        genome_summary_frame[genome] = [id_dict[i] if i in id_dict else 0 for i in genome_summary_frame.gene_id]
    return genome_summary_frame


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


def make_genome_summary(annotations, genome_summary_frame, trna_frame=None,
                        group_column='fasta', viral=False):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(fill_genome_summary_frame(annotations, group_column, genome_summary_frame.copy()))

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


##################
# Now methods for making module summary

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


def make_module_coverage_df(annotation_df, module_nets, min_cov=.5):
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
    return coverage_df.loc[coverage_df.step_coverage >= min_cov]


def make_module_coverage_summary(annotations, module_nets, min_cov=.001, group_column='fasta'):
    # go through each scaffold to check for modules
    scaffold_df_dict = dict()
    for scaffold, frame in annotations.groupby(group_column):
        scaffold_df_dict[scaffold] = make_module_coverage_df(frame, module_nets, min_cov)
    scaffold_coverage_df = pd.concat(scaffold_df_dict)
    return scaffold_coverage_df


def summarize_genomes(input_file, trna_path, output_dir, group_column, viral=False, min_cov=.001):
    # read in data
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    if trna_path is None:
        trna_frame = None
    else:
        trna_frame = pd.read_csv(trna_path, sep='\t', index_col=0)

    # get db_locs and read in dbs
    db_locs = get_database_locs()
    if 'genome_summary_form' not in db_locs or 'module_step_form' not in db_locs:
        raise ValueError('Genome_summary_frame and module_step_form must be set in order to summarize genomes.')

    # read in dbs
    genome_summary_form = pd.read_csv(db_locs['genome_summary_form'], sep='\t')
    module_steps_form = pd.read_csv(db_locs['module_step_form'], sep='\t')

    # make output folder
    mkdir(output_dir)

    # make genome summary
    genome_summary = make_genome_summary(annotations, genome_summary_form, trna_frame, group_column, viral)
    genome_summary.to_csv(path.join(output_dir, 'genome_summary.tsv'), sep='\t', index=False)

    # build module nets
    module_nets = {module: build_module_net(module_df)
                   for module, module_df in module_steps_form.groupby('module')}

    # make coverage summary
    module_coverage_summary = make_module_coverage_summary(annotations, module_nets, min_cov, group_column)
    module_coverage_summary.to_csv(path.join(output_dir, 'module_coverage_summary.tsv'), sep='\t')
