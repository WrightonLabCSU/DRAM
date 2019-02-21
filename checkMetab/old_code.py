from os import path, remove
import subprocess
import pandas as pd

from checkMetab.checkMetab import BOUTFMT6_COLUMNS


def run_mmseqs(query_db, target_db, output_db='mmseq_results.mmsdb', return_df=False, threads=10):
    subprocess.run(['mmseqs', 'search', query_db, target_db, output_db, 'tmp', '--threads', threads])
    if return_df:
        output_loc = 'output.b6'
        subprocess.run(['mmseqs', 'convertalis', query_db, target_db, output_db, output_loc, '--threads', threads])
        return pd.read_table(output_loc, header=None, names=BOUTFMT6_COLUMNS)


def get_reverse_best_hits_old(query_db, target_db, output_dir='.', bit_score_threshold=60, threads=10):
    # make query to target db
    query_target_db = path.join(output_dir, 'query_target.mmsdb')
    forward_hits = run_mmseqs(query_db, target_db, query_target_db, threads=threads)
    # get top results from target db
    forward_hits_sig = forward_hits.loc[forward_hits.bitScore > bit_score_threshold]
    forward_best_hits = dict()
    for group, frame in forward_hits_sig.groupby('qId'):
        max_bitScore = max(frame.bitScore)
        with_max = frame.loc[frame.bitScore == max_bitScore].tId
        forward_best_hits[group] = list(with_max)
    # get full best hits ids from .mmsdb_h
    short_ids_to_keep = set([id_ for ids in forward_best_hits.values() for id_ in ids])
    ids_to_keep = list()
    for line in open('%s_h' % target_db, 'rb'):
        full_id = str(line.strip())[6:-1]
        try:
            if full_id.split()[0] in short_ids_to_keep:
                ids_to_keep.append(full_id)
        except IndexError:
            print(full_id)
    ids_to_keep_loc = path.join(output_dir, 'target_ids_to_keep.txt')
    with open(ids_to_keep_loc, 'w') as f:
        f.write('%s\n' % '\n'.join(ids_to_keep))
    # filter target database for reverse search
    filtered_target_db = path.join(output_dir, 'filtered_target.mmsdb')
    subprocess.run(['mmseqs', 'createsubdb', ids_to_keep_loc, target_db, filtered_target_db])
    subprocess.run(['mmseqs', 'createsubdb', ids_to_keep_loc, '%s_h' % target_db, '%s_h' % filtered_target_db])
    subprocess.run(['cp', '%s.dbtype' % target_db, '%s.dbtype' % filtered_target_db])
    # make target to query db
    target_query_db = path.join(output_dir, 'target_query.mmsdb')
    reverse_hits = run_mmseqs(filtered_target_db, query_db, target_query_db, threads=threads)
    reverse_hits_sig = reverse_hits.loc[reverse_hits.bitScore > bit_score_threshold]
    reverse_best_hits = dict()
    for group, frame in reverse_hits_sig.groupby('qId'):
        max_bitScore = max(frame.bitScore)
        with_max = frame.loc[frame.bitScore == max_bitScore].tId
        reverse_best_hits[group] = list(with_max)
    remove(filtered_target_db)
    # analyze results
    reciprocal_best_hits = dict()
    for gene, db_ids in forward_best_hits.items():
        for db_id in db_ids:
            if db_id in reverse_best_hits:
                if gene in reverse_best_hits[db_id]:
                    reciprocal_best_hits[gene] = db_id
    return reciprocal_best_hits


def download_and_process_pfam_pfam_scan(output_dir, pfam_release='32.0', verbose=True):
    if verbose:
        print('downloading hmm file to %s' % output_dir)
    pfam_hmm_zipped = path.join(output_dir, 'Pfam-A.hmm.gz')
    subprocess.run(['wget', '-O', pfam_hmm_zipped,
                    'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.hmm.gz' % pfam_release])
    if verbose:
        print('unzipping %s' % pfam_hmm_zipped)
    subprocess.run(['gunzip', pfam_hmm_zipped])
    if verbose:
        print('pressing hmm using hmmpress')
    pfam_hmm = path.join(output_dir, 'Pfam-A.hmm')
    subprocess.run(['hmmpress', pfam_hmm])
    if verbose:
        print('downloading .hmm.dat file to %s' % output_dir)
    pfam_dat_zipped = path.join(output_dir, 'Pfam-A.hmm.dat.gz')
    subprocess.run(['wget', '-O', pfam_dat_zipped,
                    'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.hmm.dat.gz' % pfam_release])
    if verbose:
        print('unzipping %s' % pfam_dat_zipped)
    subprocess.run(['gunzip', pfam_dat_zipped])


def run_pfam_scan(gene_faa, pfam_loc, output_loc='pfam_out.txt', threads=10):
    subprocess.run(['pfam_scan.pl', '-fasta', gene_faa, '-dir', pfam_loc, '-outfile', output_loc, '-cpu', threads])



split_into_steps_space = partial(split_into_steps, split_char=' ')
split_into_steps_comma = partial(split_into_steps, split_char=',')

def parse_steps(definition, splitter=split_into_steps, flipper=False, level=''):
    steps = splitter(definition)
    ko_dict = defaultdict(dict)
    if len(steps) > 1:
        for i, step in enumerate(steps):
            if flipper:
                ko_dict.update(parse_steps(step, split_into_steps_space, False, level='%s,%s' % (level, i)))
            else:
                ko_dict.update(parse_steps(step, split_into_steps_comma, True, level='%s,%s' % (level, i)))
    else:
        level = level[1:]
        parsed_step = steps[0]
        for ko in re.split('[+|-]', parsed_step):
            ko_loc = parsed_step.find(ko)
            if ko_loc == 0:
                required = True
            elif parsed_step[ko_loc-1] == '+':
                required = True
            elif parsed_step[ko_loc-1] == '-':
                required = False
            else:
                raise ValueError('%s not expected' % parsed_step[ko_loc-1])
            ko_dict[ko] = {'step': level, 'required': required}

    return ko_dict


def generate_rows(record):
    module = record['ENTRY']
    definition = record['DEFINITION']

    pre_form = pd.DataFrame(index=['substrate', 'product', 'substrate_id', 'product_id', 'module',
                                   'path', 'step', 'gene'])
    for step, expression in enumerate(definition.split(' ')):
        # get substrate and product
        if 'COMPOUND' in record:
            substrate_id, substrate = record['COMPOUND'][step]
            product_id, product = record['COMPOUND'][step+1]
        else:
            substrate = None
            substrate_id = None
            product = None
            product_id = None
        # parse definition
        if expression.startswith('(') and expression.endswith(')'):
            expression = expression[1:-1]
        for path_number, protein in enumerate(expression.split(',')):
            path = string.ascii_uppercase[path_number]
            for ko in [j for i in protein.split('+') for j in i.split('-')]:  # TODO: Handle option subunits aka '-'
                print(module)
                protein_name = get_protein_name(ko, record['ORTHOLOGY'])
                pre_form[ko] = [substrate, product, substrate_id, product_id, module, path, step, protein_name]
    return pre_form.transpose()

def parse_steps(definition, splitter=split_into_steps, flipper=False, level=''):
    steps = splitter(definition)
    ko_dict = defaultdict(dict)
    if len(steps) > 1:
        for i, step in enumerate(steps):
            if flipper:
                ko_dict.update(parse_steps(step, split_into_steps_space, False, level='%s,%s' % (level, i)))
            else:
                ko_dict.update(parse_steps(step, split_into_steps_comma, True, level='%s,%s' % (level, i)))
    else:
        level = level[1:]
        parsed_step = steps[0]
        for ko in re.split('[+|-]', parsed_step):
            ko_loc = parsed_step.find(ko)
            if ko_loc == 0:
                required = True
            elif parsed_step[ko_loc-1] == '+':
                required = True
            elif parsed_step[ko_loc-1] == '-':
                required = False
            else:
                raise ValueError('%s not expected' % parsed_step[ko_loc-1])
            ko_dict[ko] = {'step': level, 'required': required}

    return ko_dict


def split_into_steps(definition):
    """Very fancy split on space"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char is '(':
            curr_level += 1
        if char is ')':
            curr_level -= 1
        if (curr_level == 0) and (char == ' '):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a+1:b]
        if step.startswith('(') and step.endswith(')'):
            step = step[1:-1]
        steps.append(step)
    return steps


def split_into_substeps(definition):
    """Very fancy split on comma"""
    level_counter = 0
    substep_starts = [-1]
    for i, char in enumerate(definition):
        if char == '(':
            level_counter += 1
        elif char == ')':
            level_counter -= 1
        elif char == ',' and level_counter == 0:
            substep_starts.append(i)
    substep_starts.append(len(definition))
    substeps = list()
    for a, b in pairwise(substep_starts):
        substep = definition[a+1:b]
        if substep.startswith('(') and substep.endswith(')'):
            substep = substeps[1:-1]
        substeps.append(substep)
    return substeps
