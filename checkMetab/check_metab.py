import pandas as pd
import networkx as nx
import os
from collections import defaultdict

# TODO: add RBH information to output
# TODO: add measure of redendancy of genes
# TODO: add total number of copies
# TODO: add tqdm progress bar


def build_module_net(module_df):
    # build net from a set of module paths
    num_steps = max([int(i.split(',')[0]) for i in set(module_df.path)]) + 1
    module_net = nx.DiGraph(num_steps=num_steps, module_id=list(module_df.module)[0],
                            module_name=list(module_df.module_name)[0])
    for path, frame in module_df.groupby('path'):
        split_path = [int(i) for i in path.split(',')]
        module_net.add_node(path, kos=set(frame.ko))
        # add incoming edge
        if split_path[0] == 0:
            module_net.add_edge('begin', path)
        else:
            module_net.add_edge('end_step_%s' % (split_path[0]-1), path)
        # add outgoing edge
        if split_path[0] == num_steps:
            module_net.add_edge(path, 'end')
        else:
            module_net.add_edge(path, 'end_step_%s' % split_path[0])
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


def make_coverage_df(annotation_df, module_nets, min_cov=.5):
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


def main(annotation_tsv, metabolism_db, output='.', min_cov=.5):
    # build module nets
    all_modules = pd.read_csv(metabolism_db, sep='\t')
    module_nets = {module: build_module_net(module_df) for module, module_df in all_modules.groupby('module')}
    print('Nets made')

    # go through observed kos in entire metagenome
    annotations = pd.read_csv(annotation_tsv, sep='\t', index_col=0)
    metagenome_coverage_df = make_coverage_df(annotations, module_nets, min_cov)
    metagenome_coverage_df.to_csv(os.path.join(output, 'metagenome_metab.tsv'), sep='\t')
    print('Got coverage over metagenome')

    # go through each scaffold to check for modules
    scaffold_df_dict = dict()
    for scaffold, frame in annotations.groupby('fasta'):
        scaffold_df_dict[scaffold] = make_coverage_df(frame, module_nets, min_cov)
        print('Got %s' % scaffold)
    scaffold_coverage_df = pd.concat(scaffold_df_dict)
    scaffold_coverage_df.to_csv(os.path.join(output, 'scaffold_metab.tsv'), sep='\t')
    print('Got coverage per scaffold')
