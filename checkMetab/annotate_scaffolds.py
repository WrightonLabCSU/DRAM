from skbio.io import read as read_sequence
from skbio.io import write as write_sequence
from os import path, mkdir, remove, rmdir
import subprocess
import pandas as pd
from datetime import datetime
import re
from glob import glob

# TODO: add binning information
# TODO: multiprocess prodigal by breaking up the fasta input file and then concatenate
# TODO: add ability to take into account multiple best hits as in old_code.py
# TODO: add real logging and verbose mode
# TODO: add pfam domain descriptions
# TODO: add gene locations in scaffold
# TODO: add silent mode

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd', 'tStart',
                    'tEnd', 'eVal', 'bitScore']


def download_unifref(output_dir, uniref_version='90', verbose=True):
    if verbose:
        print('downloading uniref fasta to %s' % output_dir)
    uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
    subprocess.run(['wget', '-O', uniref_fasta_zipped,
                    'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz'
                    % (uniref_version, uniref_version)], check=True)
    if verbose:
        print('unzipping %s' % uniref_fasta_zipped)
    subprocess.run(['gunzip', uniref_fasta_zipped], check=True)


def download_and_process_pfam(output_dir, pfam_release='32.0', threads=10, verbose=True):
    if verbose:
        print('downloading pfam msa to %s' % output_dir)
    pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
    subprocess.run(['wget', '-O', pfam_full_zipped,
                    'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.full.gz' % pfam_release], check=True)
    mmseq_msa = path.join(output_dir, 'pfam.mmsmsa')
    subprocess.run(['mmseqs', 'convertmsa', pfam_full_zipped, mmseq_msa], check=True)
    mmseq_profile = path.join(output_dir, 'pfam.mmspro')
    subprocess.run(['mmseqs', 'msa2profile', mmseq_msa, mmseq_profile, '--match-mode', '1', '--threads', str(threads)],
                   check=True)
    tmp_dir = path.join(output_dir, 'tmp')
    subprocess.run(['mmseqs', 'createindex', mmseq_profile, tmp_dir, '-k', '5', '-s', '7', '--threads', str(threads)],
                   check=True)


def make_mmseqs_db(fasta_loc, output_loc, create_index=False, threads=10, verbose=False):
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    subprocess.run(['mmseqs', 'createdb', fasta_loc, output_loc], check=True, stdout=stdout, stderr=stderr)
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), 'tmp')
        subprocess.run(['mmseqs', 'createindex', output_loc, tmp_dir, '--threads', str(threads)], check=True,
                       stdout=stdout, stderr=stderr)


def filter_fasta(fasta_loc, min_len=5000, output_loc=None):
    kept_seqs = (seq for seq in read_sequence(fasta_loc, format='fasta') if len(seq) > min_len)
    if output_loc is None:
        return list(kept_seqs)
    else:
        write_sequence(kept_seqs, format='fasta', into=output_loc)


def run_prodigal(fasta_loc, output_dir, verbose=False):
    output_gff = path.join(output_dir, 'genes.gff')
    output_fna = path.join(output_dir, 'genes.fna')
    output_faa = path.join(output_dir, 'genes.faa')
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    subprocess.run(['prodigal', '-i', fasta_loc, '-p', 'meta', '-f', 'gff', '-o', output_gff, '-a', output_faa, '-d',
                    output_fna], check=True, stdout=stdout, stderr=stderr)
    return output_gff, output_fna, output_faa


def get_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                  bit_score_threshold=60, threads=10, verbose=False):
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    # make query to target db
    tmp_dir = path.join(output_dir, 'tmp')
    query_target_db = path.join(output_dir, '%s_%s.mmsdb' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'search', query_db, target_db, query_target_db, tmp_dir, '--threads', str(threads)],
                   check=True, stdout=stdout, stderr=stderr)
    # filter query to target db to only best hit
    query_target_db_top = path.join(output_dir, '%s_%s.tophit.mmsdb' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'filterdb', query_target_db, query_target_db_top, '--extract-lines', '1'], check=True,
                   stdout=stdout, stderr=stderr)
    # filter query to target db to only hits with min threshold
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                    str(bit_score_threshold), '--threads', str(threads), query_target_db_top, query_target_db_top_filt],
                   check=True, stdout=stdout, stderr=stderr)
    # convert results to blast outformat 6
    forward_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'convertalis', query_db, target_db, query_target_db_top_filt, forward_output_loc,
                    '--threads', str(threads)], check=True, stdout=stdout, stderr=stderr)
    return forward_output_loc


def get_reciprocal_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                             bit_score_threshold=60, threads=10, verbose=False):
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    # make query to target db
    tmp_dir = path.join(output_dir, 'tmp')
    # create subset for second search
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    query_target_db_filt_top_swapped = path.join(output_dir, '%s_%s.minbitscore%s.tophit.swapped.mmsdb'
                                                 % (query_prefix, target_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'swapdb', query_target_db_top_filt, query_target_db_filt_top_swapped, '--threads',
                    str(threads)], check=True, stdout=stdout, stderr=stderr)
    target_db_filt = path.join(output_dir, '%s.filt.mmsdb' % target_prefix)
    subprocess.run(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, target_db, target_db_filt], check=True,
                   stdout=stdout, stderr=stderr)
    subprocess.run(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, '%s_h' % target_db,
                    '%s_h' % target_db_filt], check=True, stdout=stdout, stderr=stderr)
    # make filtered target db to query db
    target_query_db = path.join(output_dir, '%s_%s.mmsdb' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'search', target_db_filt, query_db, target_query_db, tmp_dir, '--threads', str(threads)],
                   check=True, stdout=stdout, stderr=stderr)
    # filter target to query results db
    target_query_db_filt = path.join(output_dir, '%s_%s.tophit.mmsdb' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'filterdb', target_query_db, target_query_db_filt, '--extract-lines', '1'], check=True,
                   stdout=stdout, stderr=stderr)
    # convert results to blast outformat 6
    reverse_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'convertalis', target_db_filt, query_db, target_query_db_filt, reverse_output_loc,
                    '--threads', str(threads)], check=True, stdout=stdout, stderr=stderr)
    return reverse_output_loc


def process_reciprocal_best_hits(forward_output_loc, reverse_output_loc, bit_score_threshold=350,
                                 target_prefix='target'):
    forward_hits = pd.read_csv(forward_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    forward_hits = forward_hits.set_index('qId')
    reverse_hits = pd.read_csv(reverse_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    reverse_hits = reverse_hits.set_index('qId')
    hits = pd.DataFrame(index=['%s_hit' % target_prefix, '%s_RBH' % target_prefix, '%s_identity' % target_prefix,
                               '%s_bitScore' % target_prefix, '%s_eVal' % target_prefix])
    for forward_hit, row in forward_hits.iterrows():
        rbh = False
        if row.tId in reverse_hits.index:
            if forward_hit == reverse_hits.loc[row.tId].tId and row.bitScore >= bit_score_threshold:
                rbh = True
        hits[forward_hit] = [row.tId, rbh, row.seqIdentity, row.bitScore, row.eVal]
    return hits.transpose()


def multigrep(search_terms, search_against, output='.'):  # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, 'hits.txt')
    with open(hits_file, 'w') as f:
        f.write('%s\n' % '\n'.join(search_terms))
    results = subprocess.run(['grep', '-a', '-F', '-f', hits_file, search_against], stdout=subprocess.PIPE, check=True)
    processed_results = [i.strip()[1:] for i in results.stdout.decode('ascii').split('\n')]
    remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ''}


def get_kegg_description(kegg_hits, kegg_loc):
    gene_description = list()
    ko_list = list()
    header_dict = multigrep(kegg_hits.kegg_hit, '%s_h' % kegg_loc)
    for kegg_hit in kegg_hits.kegg_hit:
        header = header_dict[kegg_hit]
        gene_description.append(header)
        kos = re.findall('(K\d\d\d\d\d)', header)
        if len(kos) == 0:
            ko_list.append('')
        else:
            ko_list.append(','.join(kos))
    new_df = pd.DataFrame([ko_list, gene_description], index=['kegg_id', 'kegg_hit'], columns=kegg_hits.index)
    return pd.concat([new_df.transpose(), kegg_hits.drop('kegg_hit', axis=1)], axis=1)


def get_uniref_description(uniref_hits, uniref_loc):
    gene_description = list()
    uniref_list = list()
    gene_taxonomy = list()
    header_dict = multigrep(uniref_hits.uniref_hit, '%s_h' % uniref_loc)
    for uniref_hit in uniref_hits.uniref_hit:
        header = header_dict[uniref_hit]
        gene_description.append(header)
        uniref_list.append(header[header.find('RepID=')+6:])
        gene_taxonomy.append(re.search('Tax=(.*?) (\S*?)=', header).group(1))
    new_df = pd.DataFrame([uniref_list, gene_description, gene_taxonomy],
                          index=['uniref_id', 'uniref_hit', 'uniref_taxonomy'],
                          columns=uniref_hits.index)
    return pd.concat([new_df.transpose(), uniref_hits.drop('uniref_hit', axis=1)], axis=1)


def run_mmseqs_pfam(query_db, pfam_profile, output_loc, output_prefix='mmpro_results', threads=10, verbose=False):
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    tmp_dir = path.join(output_loc, 'tmp')
    output_db = path.join(output_loc, '%s.mmsdb' % output_prefix)
    subprocess.run(['mmseqs', 'search', query_db, pfam_profile, output_db, tmp_dir, '-k', '5', '-s', '7', '--threads',
                    str(threads)], check=True, stdout=stdout, stderr=stderr)
    output_loc = path.join(output_loc, 'pfam_output.b6')
    subprocess.run(['mmseqs', 'convertalis', query_db, pfam_profile, output_db, output_loc], check=True, stdout=stdout,
                   stderr=stderr)
    pfam_results = pd.read_csv(output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    pfam_dict = dict()
    for gene, pfam_frame in pfam_results.groupby('qId'):
        pfam_dict[gene] = ','.join(pfam_frame.tId)
    return pd.Series(pfam_dict, name='pfam_hits')


def get_sig(tcovlen, evalue):
    if tcovlen >= 80 and evalue < 1e-5:
        return True
    elif tcovlen < 80 and evalue < 1e-3:
        return True
    else:
        return False


def run_hmmscan_dbcan(genes_faa, dbcan_loc, output_loc, verbose=False):
    """
    hmmscan --domtblout ~/dbCAN_test_1 dbCAN-HMMdb-V7.txt ~/shale_checkMetab_test/checkMetab/genes.faa
    cat ~/dbCAN_test_1 | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \
    sort -k 3,3 -k 8n -k 9n > dbCAN_test_1.good_cols.tsv
    """
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    dbcan_output = path.join(output_loc, 'dbcan_results.unprocessed.txt')
    subprocess.run(['hmmscan', '--domtblout', dbcan_output, dbcan_loc, genes_faa], check=True, stdout=stdout,
                   stderr=stderr)
    processed_dbcan_output = path.join(output_loc, 'dbcan_results.tsv')
    cmd = "cat %s | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' |" \
          "sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n > %s" % (dbcan_output, processed_dbcan_output)
    subprocess.run(cmd, shell=True, check=True)

    dbcan_res = pd.read_csv(processed_dbcan_output, sep='\t', header=None)

    columns = ['tid', 'tlen', 'qid', 'qlen', 'evalue', 'tstart', 'tend', 'qstart', 'qend']
    dbcan_res.columns = columns

    dbcan_res['tcovlen'] = dbcan_res.tend - dbcan_res.tstart

    dbcan_res['significant'] = [get_sig(row.tcovlen, row.evalue) for _, row in dbcan_res.iterrows()]

    dbcan_dict = dict()
    for gene, frame in dbcan_res[dbcan_res.significant].groupby('qid'):
        dbcan_dict[gene] = ','.join([i[:-4] for i in frame.tid])

    return pd.Series(dbcan_dict, name='cazy_hits')


def get_scaffold_and_gene(annotations):
    gene_scaffold_list = list()
    for label in annotations:
        split_label = label.split('_')
        gene_scaffold_list.append(['_'.join(split_label[:-1]), split_label[-1]])
    return pd.DataFrame(gene_scaffold_list, columns=['scaffold', 'gene_position'], index=annotations)


def get_unannotated(fasta_loc, annotations):
    return [seq.metadata['id'] for seq in read_sequence(fasta_loc, format='fasta')
            if seq.metadata['id'] not in annotations]


def assign_grades(annotations):
    grades = dict()
    for gene, row in annotations.iterrows():
        if row.kegg_RBH is True:
            grade = 'A'
        elif row.uniref_RBH is True:
            grade = 'B'
        elif not pd.isna(row.kegg_hit) or not pd.isna(row.uniref_hit):
            grade = 'C'
        elif not pd.isna(row.pfam_hits):
            grade = 'D'
        else:
            grade = 'E'
        grades[gene] = grade
    return pd.Series(grades, name='grade')


def generate_annotated_fasta(input_fasta, annotations, verbosity='short', name=None):
    """verbosity should be short or long"""
    for seq in read_sequence(input_fasta, format='fasta'):
        annotation = annotations.loc[seq.metadata['id']]
        annotation_str = 'grade: %s' % annotation.grade
        if verbosity == 'short':
            if (annotation.grade == 'A') or (annotation.grade == 'C' and not pd.isna(annotation.kegg_hit)):
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'kegg')
            if annotation.grade == 'B' or (annotation.grade == 'C' and not pd.isna(annotation.uniref_hit)):
                annotation_str += '; %s (db=%s)' % (annotation.uniref_hit, 'uniref')
            if annotation.grade == 'D':
                annotation_str += '; %s (db=%s)' % (annotation.pfam_hits, 'pfam')
        elif verbosity == 'long':
            if not pd.isna(annotation.kegg_hit):
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'kegg')
            if not pd.isna(annotation.uniref_hit):
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'uniref')
            if not pd.isna(annotation.pfam_hits):
                annotation_str += '; %s (db=%s)' % (annotation.pfam_hits, 'pfam')
        else:
            raise ValueError('%s is not a valid verbosity level for annotation summarization' % verbosity)
        if name is not None:
            seq.metadata['id'] = '%s_%s' % (name, seq.metadata['id'])
        seq.metadata['description'] = annotation_str
        yield seq


def create_annotated_fasta(input_fasta, annotations, output_fasta, verbosity='short', name=None):
    write_sequence(generate_annotated_fasta(input_fasta, annotations, verbosity, name),
                   format='fasta', into=output_fasta)


def generate_renamed_fasta(input_fasta, prefix):
    for seq in read_sequence(input_fasta, format='fasta'):
        seq.metadata['id'] = '%s_%s' % (prefix, seq.metadata['id'])
        yield seq


def rename_fasta(input_fasta, output_fasta, prefix):
    write_sequence(generate_renamed_fasta(input_fasta, prefix), format='fasta', into=output_fasta)


def rename_gff(input_gff, output_gff, prefix):
    with open(input_gff) as f:
        with open(output_gff, 'w') as o:
            for line in f:
                if not line.startswith('#') and not line.startswith('\n'):
                    old_scaffold = line.strip().split('\t')[0]
                    line = '%s_%s' % (prefix, line)
                    match = re.search('ID=\d*_\d*;', line)
                    gene_number = match.group().split('_')[-1][:-1]
                    line = re.sub('ID=\d*_\d*;', 'ID=%s_%s_%s;' % (prefix, old_scaffold, gene_number), line)
                o.write(line)


def merge_files(files_to_merge, outfile):
    with open(outfile, 'w') as outfile_handle:
        for file in glob(files_to_merge):
            with open(file) as f:
                outfile_handle.write(f.read())


def merge_gtfs(gtf_files, outfile):
    gtf_files = glob(gtf_files)
    with open(outfile, 'w') as f:
        f.write(open(gtf_files[0]).readline())
        for gtf in gtf_files:
            content = ''.join(open(gtf).readlines()[1:])
            f.write(content)


def main(fasta_glob_str, kegg_loc, uniref_loc, pfam_loc, dbcan_loc, output_dir='.', min_size=5000,
         bit_score_threshold=60, rbh_bit_score_threshold=350, keep_tmp_dir=True, threads=10, verbose=True):
    # set up
    start_time = datetime.now()
    fasta_locs = glob(fasta_glob_str)
    if len(fasta_locs) == 0:
        raise ValueError('Given fasta locations returns no paths: %s')
    else:
        print('%s fasta found' % len(fasta_locs))
    if not path.isfile(kegg_loc):
        raise ValueError('KEGG mmsdb does not exist: %s' % kegg_loc)
    if not path.isfile(uniref_loc):
        raise ValueError('UniRef mmsdb does not exist: %s' % uniref_loc)
    if not path.isfile(uniref_loc):
        raise ValueError('PFam mmspro does not exist: %s' % pfam_loc)
    mkdir(output_dir)
    tmp_dir = path.join(output_dir, 'working_dir')
    mkdir(tmp_dir)

    annotations_list = list()
    for fasta_loc in fasta_locs:
        fasta_name = path.splitext(path.basename(fasta_loc))[0]
        print('%s: Annotating %s' % (str(datetime.now()-start_time), fasta_name))
        fasta_dir = path.join(tmp_dir, fasta_name)
        mkdir(fasta_dir)

        # first step filter fasta
        print('%s: Filtering fasta' % str(datetime.now()-start_time))
        filtered_fasta = path.join(fasta_dir, 'filtered_fasta.fa')
        filter_fasta(fasta_loc, min_size, filtered_fasta)

        # call genes with prodigal
        print('%s: Calling genes with prodigal' % str(datetime.now()-start_time))
        gene_gff, gene_fna, gene_faa = run_prodigal(filtered_fasta, fasta_dir, verbose=verbose)

        # run reciprocal best hits from kegg and uniref
        print('%s: Turning genes from prodigal to mmseqs2 db' % str(datetime.now()-start_time))
        query_db = path.join(fasta_dir, 'gene.mmsdb')
        make_mmseqs_db(gene_faa, query_db, create_index=True, threads=threads, verbose=verbose)

        annotation_list = list()

        if kegg_loc is not None:
            print('%s: Getting forward best hits from KEGG' % str(datetime.now()-start_time))
            forward_kegg_hits = get_best_hits(query_db, kegg_loc, fasta_dir, 'gene', 'kegg', bit_score_threshold,
                                              threads, verbose=verbose)
            print('%s: Getting reverse best hits from KEGG' % str(datetime.now()-start_time))
            reverse_kegg_hits = get_reciprocal_best_hits(query_db, kegg_loc, fasta_dir, 'gene', 'kegg',
                                                         bit_score_threshold, threads, verbose=verbose)
            kegg_hits = process_reciprocal_best_hits(forward_kegg_hits, reverse_kegg_hits, rbh_bit_score_threshold,
                                                     'kegg')
            kegg_hits = get_kegg_description(kegg_hits, kegg_loc)
            annotation_list.append(kegg_hits)

        if uniref_loc is not None:
            print('%s: Getting forward best hits from UniRef' % str(datetime.now()-start_time))
            forward_uniref_hits = get_best_hits(query_db, uniref_loc, fasta_dir, 'gene', 'uniref', bit_score_threshold,
                                                threads, verbose=verbose)
            print('%s: Getting reverse best hits from UniRef' % str(datetime.now()-start_time))
            reverse_uniref_hits = get_reciprocal_best_hits(query_db, uniref_loc, fasta_dir, 'gene', 'uniref',
                                                           bit_score_threshold, threads, verbose=verbose)
            uniref_hits = process_reciprocal_best_hits(forward_uniref_hits, reverse_uniref_hits,
                                                       rbh_bit_score_threshold, 'uniref')
            uniref_hits = get_uniref_description(uniref_hits, uniref_loc)
            annotation_list.append(uniref_hits)

        # run pfam scan
        if pfam_loc is not None:
            print('%s: Getting hits from pfam' % str(datetime.now()-start_time))
            pfam_hits = run_mmseqs_pfam(query_db, pfam_loc, fasta_dir, output_prefix='pfam', threads=threads,
                                        verbose=verbose)
            annotation_list.append(pfam_hits)

        # use hmmer to detect cazy ids using dbCAN
        if dbcan_loc is not None:
            print('%s: Getting hits from dbCAN' % str(datetime.now()-start_time))
            dbcan_hits = run_hmmscan_dbcan(gene_faa, dbcan_loc, fasta_dir, verbose=verbose)
            annotation_list.append(dbcan_hits)

        # merge dataframes
        print('%s: Finishing up results' % str(datetime.now()-start_time))
        annotations = pd.concat(annotation_list, axis=1, sort=False)

        # get scaffold data and assign grades
        if uniref_loc is not None and kegg_loc is not None:
            grades = assign_grades(annotations)
            annotations = pd.concat([grades, annotations], axis=1)
        annotations = pd.concat([get_scaffold_and_gene(annotations.index), annotations], axis=1)

        # add unknowns
        unannotated_genes = get_unannotated(gene_faa, annotations.index)
        unannotated = get_scaffold_and_gene(unannotated_genes)
        unannotated['grade'] = 'E'
        annotations = pd.concat([unannotated, annotations], sort=False)

        # generate fna and faa output files with uniref annotations
        annotated_fna = path.join(fasta_dir, 'genes.annotated.fna')
        create_annotated_fasta(gene_fna, annotations, annotated_fna, name=fasta_name)
        annotated_faa = path.join(fasta_dir, 'genes.annotated.faa')
        create_annotated_fasta(gene_faa, annotations, annotated_faa, name=fasta_name)
        renamed_scaffolds = path.join(fasta_dir, 'scaffolds.annotated.fa')
        rename_fasta(filtered_fasta, renamed_scaffolds, prefix=fasta_name)
        renamed_gffs = path.join(fasta_dir, 'genes.annotated.gff')
        rename_gff(gene_gff, renamed_gffs, prefix=fasta_name)

        # add fasta name to frame and index, append to list
        annotations.insert(0, 'fasta', fasta_name)
        annotations.index = annotations.fasta + '_' + annotations.index
        annotations_list.append(annotations)

    # merge annotation dicts
    all_annotations = pd.concat(annotations_list, sort=False)
    all_annotations = all_annotations.sort_values(['fasta', 'scaffold', 'gene_position'])
    all_annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')

    # merge gene files
    merge_files(path.join(tmp_dir, '*', '*.annotated.fna'), path.join(output_dir, 'genes.fna'))
    merge_files(path.join(tmp_dir, '*', '*.annotated.faa'), path.join(output_dir, 'genes.faa'))
    merge_files(path.join(tmp_dir, '*', 'scaffolds.annotated.fa'), path.join(output_dir, 'scaffolds.fna'))
    merge_gtfs(path.join(tmp_dir, '*', 'genes.annotated.gff'), path.join(output_dir, 'genes.gff'))

    # clean up
    if not keep_tmp_dir:
        rmdir(tmp_dir)

    print("%s: Completed" % str(datetime.now()-start_time))
