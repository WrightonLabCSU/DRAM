from skbio.io import read as read_sequence
from skbio.io import write as write_sequence
from os import path, mkdir
import subprocess
import pandas as pd
from datetime import datetime
import re

# TODO: add binning information
# TODO: multiprocess prodigal by breaking up the fasta input file and then concatenate
# TODO: add ability to take into account multiple best hits as in old_code.py
# TODO: add real logging and verbose mode

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd', 'tStart',
                    'tEnd', 'eVal', 'bitScore']


def download_unifref(output_dir, uniref_version='90', verbose=True):
    if verbose:
        print('downloading uniref fasta to %s' % output_dir)
    uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
    subprocess.run(['wget', '-O', uniref_fasta_zipped,
                    'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz'
                    % (uniref_version, uniref_version)])
    if verbose:
        print('unzipping %s' % uniref_fasta_zipped)
    subprocess.run(['gunzip', uniref_fasta_zipped])


def download_and_process_pfam(output_dir, pfam_release='32.0', threads=10, verbose=True):
    if verbose:
        print('downloading pfam msa to %s' % output_dir)
    pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
    subprocess.run(['wget', '-O', pfam_full_zipped,
                    'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.full.gz' % pfam_release])
    mmseq_msa = path.join(output_dir, 'pfam.mmsmsa')
    subprocess.run(['mmseqs', 'convertmsa', pfam_full_zipped, mmseq_msa])
    mmseq_profile = path.join(output_dir, 'pfam.mmspro')
    subprocess.run(['mmseqs', 'msa2profile', mmseq_msa, mmseq_profile, '--match-mode', '1', '--threads', str(threads)])
    tmp_dir = path.join(output_dir, 'tmp')
    subprocess.run(['mmseqs', 'createindex', mmseq_profile, tmp_dir, '-k', '5', '-s', '7', '--threads', str(threads)])


def make_mmseqs_db(fasta_loc, output_loc, create_index=False, threads=10):
    subprocess.run(['mmseqs', 'createdb', fasta_loc, output_loc])
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), 'tmp')
        subprocess.run(['mmseqs', 'createindex', output_loc, tmp_dir, '--threads', str(threads)])


def filter_fasta(fasta_loc, min_len=5000, output_loc=None):
    kept_seqs = (seq for seq in read_sequence(fasta_loc, format='fasta') if len(seq) > min_len)
    if output_loc is None:
        return kept_seqs
    else:
        write_sequence(kept_seqs, format='fasta', into=output_loc)


def run_prodigal(fasta_loc, output_dir):
    output_gff = path.join(output_dir, 'genes.gff')
    output_fna = path.join(output_dir, 'genes.fna')
    output_faa = path.join(output_dir, 'genes.faa')
    subprocess.run(['prodigal', '-i', fasta_loc, '-p', 'meta', '-f', 'gff',
                    '-o', output_gff, '-a', output_faa, '-d', output_fna])
    return output_gff, output_fna, output_faa


def get_reciprocal_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                             bit_score_threshold=60, threads=10):
    # make query to target db
    tmp_dir = path.join(output_dir, 'tmp')
    query_target_db = path.join(output_dir, '%s_%s.mmsdb' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'search', query_db, target_db, query_target_db, tmp_dir, '--threads', str(threads)])
    # filter query to target db to only best hit
    query_target_db_top = path.join(output_dir, '%s_%s.tophit.mmsdb' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'filterdb', query_target_db, query_target_db_top, '--extract-lines', '1'])
    # filter query to target db to only hits with min threshold
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%smmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                    str(bit_score_threshold), '--threads', str(threads), query_target_db_top, query_target_db_top_filt])
    # create subset for second search
    query_target_db_filt_top_swapped = path.join(output_dir, '%s_%s.minbitscore%s.tophit.swapped.mmsdb'
                                                 % (query_prefix, target_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'swapdb', query_target_db_top_filt, query_target_db_filt_top_swapped, '--threads',
                    str(threads)])
    target_db_filt = path.join(output_dir, '%s.filt.mmsdb' % target_prefix)
    subprocess.run(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, target_db, target_db_filt])
    subprocess.run(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, '%s_h' % target_db,
                    '%s_h' % target_db_filt])
    # make filtered target db to query db
    target_query_db = path.join(output_dir, '%s_%s.mmsdb' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'search', target_db_filt, query_db, target_query_db, tmp_dir, '--threads', str(threads)])
    # filter target to query results db
    target_query_db_filt = path.join(output_dir, '%s_%s.tophit.mmsdb' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'filterdb', target_query_db, target_query_db_filt, '--extract-lines', '1'])
    # convert results to blast outformat 6
    forward_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'convertalis', query_db, target_db, query_target_db_top_filt, forward_output_loc,
                    '--threads', str(threads)])
    reverse_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'convertalis', target_db_filt, query_db, target_query_db_filt, reverse_output_loc,
                    '--threads', str(threads)])
    return forward_output_loc, reverse_output_loc


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


def get_kegg_description(kegg_hits, kegg_loc):
    gene_description = list()
    ko_list = list()
    for kegg_hit in kegg_hits.kegg_hit:
        result = subprocess.run(['grep', '-a', kegg_hit, '%s_h' % kegg_loc], capture_output=True)
        header = result.stdout.decode('ascii').strip()[1:]
        gene_description.append(header)
        kos = re.findall('(K\d\d\d\d\d)', header)
        if len(kos) == 0:
            ko_list.append('')
        else:
            ko_list.append(','.join(kos))
    kegg_hits['kegg_hit'] = gene_description
    kegg_hits['kegg_ko'] = ko_list
    return kegg_hits


def get_uniref_description(uniref_hits, uniref_loc):
    gene_description = list()
    uniref_list = list()
    for uniref_hit in uniref_hits.uniref_hit:
        result = subprocess.run(['grep', '-a', uniref_hit, '%s_h' % uniref_loc], capture_output=True)
        header = result.stdout.decode('ascii').strip()[1:]
        gene_description.append(header)
        uniref_list.append(header[header.find('RepID=')+6:])
    uniref_hits['uniref_hit'] = gene_description
    uniref_hits['uniref_id'] = uniref_list
    return uniref_hits


def run_mmseqs_pfam(query_db, pfam_profile, output_loc, output_prefix='mmpro_results', threads=10):
    tmp_dir = path.join(output_loc, 'tmp')
    output_db = path.join(output_loc, '%s.mmsdb' % output_prefix)
    subprocess.run(['mmseqs', 'search', query_db, pfam_profile, output_db, tmp_dir, '-k', '5', '-s', '7', '--threads',
                    str(threads)])
    output_loc = path.join(output_loc, 'pfam_output.b6')
    subprocess.run(['mmseqs', 'convertalis', query_db, pfam_profile, output_db, output_loc])
    pfam_results = pd.read_csv(output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    pfam_dict = dict()
    for gene, pfam_frame in pfam_results.groupby('qId'):
        pfam_dict[gene] = ','.join(pfam_frame.tId)
    return pd.Series(pfam_dict, name='pfam_hits')


def get_scaffold(annotations):
    return pd.Series(['_'.join(label.split('_')[:-1]) for label in annotations.index], index=annotations.index,
                     name='scaffold')


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


def main(fasta_loc, kegg_loc, uniref_loc, pfam_loc, output_dir='.', min_size=5000, bit_score_threshold=60,
         rbh_bit_score_threshold=350, threads=10):
    start_time = datetime.now()
    mkdir(output_dir)
    # first step filter fasta
    print('Filtering fasta')
    filtered_fasta = path.join(output_dir, 'filtered_fasta.fa')
    filter_fasta(fasta_loc, min_size, filtered_fasta)
    # call genes with prodigal
    print('Calling genes with prodigal')
    gene_gff, gene_fna, gene_faa = run_prodigal(filtered_fasta, output_dir)
    # run reverse best hits from kegg and uniref
    print('Turning genes from prodigcal to mmseqs2 db')
    query_db = path.join(output_dir, 'gene.mmsdb')
    make_mmseqs_db(gene_faa, query_db, create_index=True, threads=threads)
    print('Getting reverse best hits from KEGG')
    forward_kegg_hits, reverse_kegg_hits = get_reciprocal_best_hits(query_db, kegg_loc, output_dir, 'gene', 'kegg',
                                                                    bit_score_threshold, threads)
    kegg_hits = process_reciprocal_best_hits(forward_kegg_hits, reverse_kegg_hits, rbh_bit_score_threshold, 'kegg')
    kegg_hits = get_kegg_description(kegg_hits, kegg_loc)
    print('Getting reverse best hits from UniRef')
    forward_uniref_hits, reverse_uniref_hits = get_reciprocal_best_hits(query_db, uniref_loc, output_dir, 'gene',
                                                                        'uniref', bit_score_threshold, threads)
    uniref_hits = process_reciprocal_best_hits(forward_uniref_hits, reverse_uniref_hits, rbh_bit_score_threshold,
                                               'uniref')
    uniref_hits = get_uniref_description(uniref_hits, uniref_loc)
    # run pfam scan
    print('Getting hits from pfam')
    pfam_hits = run_mmseqs_pfam(query_db, pfam_loc, output_dir, output_prefix='pfam', threads=threads)
    # merge dataframes
    print('Finishing up results')
    annotations = pd.concat([kegg_hits, uniref_hits, pfam_hits], axis=1)
    # get scaffold data, assign grade and output
    annotations = pd.concat([get_scaffold(annotations), assign_grades(annotations), annotations], axis=1)
    annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t', index_label='gene')
    print("Runtime: %s" % str(datetime.now()-start_time))
