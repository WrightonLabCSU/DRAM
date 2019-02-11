from skbio.io import read as read_sequence
from skbio.io import write as write_sequence
from os import path, mkdir
import subprocess
import pandas as pd

# TODO: add binning information
# TODO: multiprocess prodigal by breaking up the fasta input file and then concatenate
# TODO: add ability to take into account multiple best hits as in old code

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
    subprocess.run(['mmseqs', 'createindex', mmseq_profile, 'tmp', '-k', '5', '-s', '7', '--threads', str(threads)])


def make_mmseqs_db(fasta_loc, output_loc, create_index=False, threads=10):
    subprocess.run(['mmseqs', 'createdb', fasta_loc, output_loc])
    if create_index:
        subprocess.run(['mmseqs', 'createindex', output_loc, 'tmp', '--threads', str(threads)])


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


def get_reverse_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                          bit_score_threshold=60, threads=10):
    # make query to target db
    query_target_db = path.join(output_dir, '%s_%s.mmsdb' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'search', query_db, target_db, query_target_db, 'tmp', '--threads', str(threads)])
    # filter query to target db to only best hit
    query_target_db_top = path.join(output_dir, '%s_%s.tophit.mmsdb' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'filterdb', query_target_db, query_target_db_top, '--extract-lines', 1])
    # filter query to target db to only hits with min threshold
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%smmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                    bit_score_threshold, query_target_db_top, query_target_db_top_filt])
    # create subset for second search
    query_target_db_filt_top_swapped = path.join(output_dir, '%s_%s.minbitscore%s.tophit.swapped.mmsdb'
                                                 % (query_prefix, target_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'swapdb', query_target_db_top_filt, query_target_db_filt_top_swapped])
    target_db_filt = path.join(output_dir, '%s.filt.mmsdb' % target_prefix)
    subprocess.run(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, target_db, target_db_filt])
    subprocess.run(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, '%s_h' % target_db,
                    '%s_h' % target_db_filt])
    # make filtered target db to query db
    target_query_db = path.join(output_dir, '%s_%s.mmsdb' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs,' 'search', target_db_filt, query_db, target_query_db, 'tmp', '--threads', str(threads)])
    # filter target to query results db
    target_query_db_filt = path.join(output_dir, '%s_%s.tophit.mmsdb' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'filterdb', target_query_db, target_query_db_filt, '--extract-lines', 1])
    # get results
    forward_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (query_prefix, target_prefix))
    subprocess.run(['mmseqs', 'convertalis', query_db, target_db, query_target_db_top_filt, forward_output_loc,
                    '--threads', str(threads)])
    forward_hits = pd.read_table(forward_output_loc, header=None, names=BOUTFMT6_COLUMNS)
    forward_hits = forward_hits.set_index('qId')
    reverse_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (target_prefix, query_prefix))
    subprocess.run(['mmseqs', 'convertalis', target_db_filt, query_db, target_query_db_filt, reverse_output_loc,
                    '--threads', str(threads)])
    reverse_hits = pd.read_table(reverse_output_loc, header=None, names=BOUTFMT6_COLUMNS)
    reverse_hits = reverse_hits.set_index('qId')
    hits = pd.read_table(index=['%s_hit' % target_prefix, '%s_RBH' % target_prefix])
    for forward_hit, row in forward_hits.iterrows():
        rbh = False
        if row.tId in reverse_hits.index:
            if forward_hit == reverse_hits.loc[row.tId].tId:
                rbh = True
        hits[forward_hit] = [row.tId, rbh]
    return hits


def run_mmseqs_pfam(query_db, pfam_profile, output_loc, output_prefix='mmpro_results', bit_score_threshold=60,
                    threads=10):
    output_db = path.join(output_loc, '%s.mmsdb' % output_prefix)
    subprocess.run(['mmseqs', 'search', query_db, pfam_profile, output_db, 'tmp', '-k', '5', '-s', '7', '--threads',
                    str(threads)])
    # filter to remove bad hits
    output_db_filt = path.join(output_loc, '%s.minbitscore%s.mmsdb' % (output_prefix, bit_score_threshold))
    subprocess.run(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                    bit_score_threshold, output_db, output_db_filt])
    output_loc = 'pfam_output.b6'
    subprocess.run(['mmseqs', 'convertalis', query_db, pfam_profile, output_db_filt, output_loc])
    pfam_results = pd.read_table(output_loc, header=None, names=BOUTFMT6_COLUMNS)
    pfam_dict = dict()
    for gene, pfam_frame in pfam_results.groupby('qId'):
        pfam_dict[gene] = ','.join(pfam_frame.tId)
    return pd.Series(pfam_dict)


def assign_grades(annotations):
    grades = list()
    for gene, row in annotations.iterrows():
        if row.kegg_RBH:
            grade = 'A'
        elif row.uniref_RBH:
            grade = 'B'
        elif row.kegg is not None:
            grade = 'C'
        elif row.uniref is not None:
            grade = 'D'
        elif row.pfam_hits is not None:
            grade = 'E'
        else:
            grade = 'F'
        grades.append(grade)
    annotations['grade'] = grades
    return annotations


def main(fasta_loc, pfam_loc, uniref_loc, kegg_loc, output_dir='.', min_size=5000, bit_score_threshold=60, threads=10):
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
    make_mmseqs_db(filtered_fasta, query_db, create_index=True, threads=threads)
    print('Getting reverse best hits from KEGG')
    kegg_hits = get_reverse_best_hits(filtered_fasta, kegg_loc, output_dir, 'gene', 'kegg', bit_score_threshold,
                                      threads)
    print('Getting reverse best hits from UniRef')
    uniref_hits = get_reverse_best_hits(filtered_fasta, uniref_loc, output_dir, 'gene', 'uniref', bit_score_threshold,
                                        threads)
    # run pfam scan
    print('Getting hits from pfam')
    pfam_hits = run_mmseqs_pfam(gene_faa, pfam_loc, output_dir, output_prefix='pfam', bit_score_threshold=60,
                                threads=threads)
    # merge dataframes
    print('Finishing up results')
    annotations = pd.concat([kegg_hits, uniref_hits])
    annotations['pfam_hits'] = pfam_hits
    # assign grade and output
    annotations = assign_grades(annotations)
    annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')
