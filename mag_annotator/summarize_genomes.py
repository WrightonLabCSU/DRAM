import pandas as pd
from collections import Counter, defaultdict
import networkx as nx
from os import path, mkdir

from mag_annotator.utils import get_database_locs

# TODO: add RBH information to output
# TODO: add measure of redendancy of genes
# TODO: add total number of copies
# TODO: add tqdm progress bar


def summarize_kos_from_summary_frame(annotations, group_column, genome_summary_frame):
    grouped = annotations.groupby(group_column)
    for genome, frame in grouped:
        kos = Counter([j for i in frame.kegg_id if type(i) is str for j in i.split(',')])
        genome_summary_frame[genome] = [kos[i] if i in kos else 0 for i in genome_summary_frame.KO]
    return genome_summary_frame


def summarize_cazys(annotations, groupby_column, genome_summary_frame_columns):
    cazy_classes = ['AA', 'CBM', 'CE', 'GH', 'GT', 'PL', 'SLH', 'cohesin', 'dockerin']
    cazy_lists = [i for i in list(annotations.cazy_hits) if type(i) is not float]
    cazy_set = sorted({j for i in cazy_lists for j in i.split(',')})
    rows = list()
    for cazy_class in cazy_classes:
        for cazy in cazy_set:
            if cazy_class in cazy:
                rows.append(('na', cazy, '%s ?' % cazy_class, 'na', 'na', 'CAZY'))
    cazy_summary_frame = pd.DataFrame(rows, columns=genome_summary_frame_columns)

    # get cazy abundances per genome
    grouped_cazy = annotations.groupby(groupby_column)
    for genome, frame in grouped_cazy:
        kos = Counter([j for i in frame.cazy_hits if type(i) is str for j in i.split(',')])
        cazy_summary_frame[genome] = [kos[i] if i in kos else 0 for i in cazy_summary_frame.KO]
    return cazy_summary_frame


def summarize_trnas(trnas_frame, groupby_column='fasta'):
    df_dict = dict()
    for virus, frame in trnas_frame.groupby(groupby_column):
        df_dict[virus] = Counter(frame.Type)
    trnas_summary_frame = pd.DataFrame.from_dict(df_dict)
    trnas_summary_frame = trnas_summary_frame.fillna(0)
    return trnas_summary_frame


def make_genome_summary(annotations, genome_summary_frame, trna_frame=None,
                        group_column='fasta', viral=False):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(summarize_kos_from_summary_frame(annotations, group_column, genome_summary_frame.copy()))

    # add cazys
    summary_frames.append(summarize_cazys(annotations, group_column, genome_summary_frame.columns))

    # add tRNAs
    if trna_frame is not None:
        summary_frames.append(summarize_trnas(trna_frame, group_column))

    # merge summary frames
    summarized_genomes = pd.concat(summary_frames)

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
    genome_summary.to_csv(path.join(output_dir, 'genome_summary.tsv'), sep='\t')

    # build module nets
    module_nets = {module: build_module_net(module_df)
                   for module, module_df in module_steps_form.groupby('module')}

    # make coverage summary
    module_coverage_summary = make_module_coverage_summary(annotations, module_nets, min_cov, group_column)
    module_coverage_summary.to_csv(path.join(output_dir, 'module_coverage_summary.tsv'), sep='\t')
