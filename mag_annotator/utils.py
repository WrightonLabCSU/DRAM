import re
import subprocess
from os import path, stat
from urllib.request import urlopen, urlretrieve
from urllib.error import URLError
import pandas as pd
import logging
from typing import Callable

HMMSCAN_ALL_COLUMNS = ['query_id', 'query_ascession', 'query_length', 'target_id', 'target_ascession', 'target_length',
                       'full_evalue', 'full_score', 'full_bias', 'domain_number', 'domain_count', 'domain_cevalue',
                       'domain_ievalue', 'domain_score', 'domain_bias', 'target_start', 'target_end', 'alignment_start',
                       'alignment_end', 'query_start', 'query_end', 'accuracy', 'description']
HMMSCAN_COLUMN_TYPES = [str, str, int, str, str, int, float, float, float, int, int, float, float, float, float, int,
                        int, int, int, int, int, float, str]
BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd', 'tStart',
                    'tEnd', 'eVal', 'bitScore']


def download_file(url: str, output_file: str, logger: logging.Logger, alt_urls: list = None, verbose = True):
    # TODO: catching error 4 and give error message to retry or retry automatically
    links = [url] if alt_urls is None else [url] + alt_urls
    for l in links: 
        if verbose:
            print('downloading %s' % url)
        try:
            urlretrieve(l, output_file)
            return
        except BaseException as error:
            # BaseException is good http was to exact
            logger.warning(f"Something went wrong with the download of the url: {l}")
            logger.warning(error)
    raise URLError("DRAM whas not able to download a key database, check the logg for details")
    # run_process(['wget', '-O', output_file, url], verbose=verbose)


def setup_logger(logger, *log_file_paths, level=logging.INFO):
    logger.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    for log_file_path in log_file_paths:
        fh = logging.FileHandler(log_file_path)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)


def run_process(command, logger, shell:bool=False, capture_stdout:bool=True, save_output:str=None, 
                check:bool=False, stop_on_error:bool=True, verbose:bool=False) -> str:
    """
    Standardization of parameters for using subprocess.run, provides verbose mode and option to run via shell
    """
    # TODO just remove check
    try:
        results = subprocess.run(command, check=check, shell=shell,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as error:
        logger.critical(f'The subcommand {command} experienced an error')
        if stop_on_error:
            raise error
    if results.returncode != 0:
        logger.critical(f'The subcommand {command} experienced an error: {results.stderr}')
        logging.debug(results.stdout)
        if stop_on_error:
           raise subprocess.SubprocessError(f"The subcommand {' '.join(command)} experienced an error, see the log for more info.")
    if save_output is not None:
        with open(save_output, 'w') as out:
            out.write(results.stdout)
    if capture_stdout:
        return results.stdout


# TODO: refactor following to methods to a shared run hmm step and individual get description steps
def parse_hmmsearch_domtblout(file):
    df_lines = list()
    for line in open(file):
        if not line.startswith('#'):
            line = line.split()
            line = line[:22] + [' '.join(line[22:])]
            df_lines.append(line)
    hmmsearch_frame = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    for i, column in enumerate(hmmsearch_frame.columns):
        hmmsearch_frame[column] = hmmsearch_frame[column].astype(HMMSCAN_COLUMN_TYPES[i])
    return hmmsearch_frame


def make_mmseqs_db(fasta_loc, output_loc, logger, create_index=True, threads=10, verbose=False):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2"""
    run_process(['mmseqs', 'createdb', fasta_loc, output_loc], logger, verbose=verbose)
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), 'tmp')
        run_process(['mmseqs', 'createindex', output_loc, tmp_dir, '--threads', str(threads)], logger, verbose=verbose)


def multigrep(search_terms, search_against, logger, split_char='\n', output='.'):
    # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, 'hits.txt')
    with open(hits_file, 'w') as f:
        f.write('%s\n' % '\n'.join(search_terms))
    results = run_process(['grep', '-a', '-F', '-f', hits_file, search_against], logger, capture_stdout=True, verbose=False)
    processed_results = [i.strip() for i in results.strip().split(split_char)
                         if len(i) > 0]
    # remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ''}


def merge_files(files_to_merge, outfile, has_header=False):
    """It's in the name, if has_header assumes all files have the same header"""
    with open(outfile, 'w') as outfile_handle:
        if has_header:
            outfile_handle.write(open(files_to_merge[0]).readline())
        for file in files_to_merge:
            with open(file) as f:
                if has_header:
                    _ = f.readline()
                outfile_handle.write(f.read())


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text  # or whatever


def remove_suffix(text, suffix):
    if text.endswith(suffix):
        return text[:-1*len(suffix)]
    return text  # or whatever


def get_ordered_uniques(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x) or pd.isna(x))]


def get_best_hits(query_db, target_db, logger, output_dir='.', query_prefix='query', target_prefix='target',
                  bit_score_threshold=60, threads=10, verbose=False):
    """Uses mmseqs2 to do a blast style search of a query db against a target db, filters to only include best hits
    Returns a file location of a blast out format 6 file with search results
    """
    # make query to target db
    tmp_dir = path.join(output_dir, 'tmp')
    query_target_db = path.join(output_dir, '%s_%s.mmsdb' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'search', query_db, target_db, query_target_db, tmp_dir, '--threads', str(threads)],
                 logger, verbose=verbose)
    # filter query to target db to only best hit
    query_target_db_top = path.join(output_dir, '%s_%s.tophit.mmsdb' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'filterdb', query_target_db, query_target_db_top, '--extract-lines', '1'], logger,
                verbose=verbose)
    # filter query to target db to only hits with min threshold
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    run_process(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                 str(bit_score_threshold), '--threads', str(threads), query_target_db_top, query_target_db_top_filt],
                logger, verbose=verbose)
    # convert results to blast outformat 6
    forward_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'convertalis', query_db, target_db, query_target_db_top_filt, forward_output_loc,
                 '--threads', str(threads)], logger, verbose=verbose)
    return forward_output_loc


def get_reciprocal_best_hits(query_db, target_db, logger, output_dir='.', query_prefix='query', target_prefix='target',
                             bit_score_threshold=60, rbh_bit_score_threshold=350, threads=10, verbose=False):
    """Take results from best hits and use for a reciprocal best hits search"""
    # TODO: Make it take query_target_db as a parameter
    # create subset for second search
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))  # I DON'T LIKE THIS
    query_target_db_filt_top_swapped = path.join(output_dir, '%s_%s.minbitscore%s.tophit.swapped.mmsdb'
                                                 % (query_prefix, target_prefix, bit_score_threshold))
    # swap queries and targets in results database
    run_process(['mmseqs', 'swapdb', query_target_db_top_filt, query_target_db_filt_top_swapped, '--threads',
                 str(threads)], logger, verbose=verbose)
    target_db_filt = path.join(output_dir, '%s.filt.mmsdb' % target_prefix)
    # create a subdatabase of the target database with the best hits as well as the index of the target database
    run_process(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, target_db, target_db_filt],
                logger, verbose=verbose)
    run_process(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, '%s_h' % target_db,
                 '%s_h' % target_db_filt], logger, verbose=verbose)

    return get_best_hits(target_db_filt, query_db, logger, output_dir, target_prefix, query_prefix, rbh_bit_score_threshold,
                         threads, verbose)


def run_hmmscan(genes_faa:str, db_loc:str, db_name:str, output_loc:str, formater:Callable,
                logger:logging.Logger, threads:int=2, db_handler=None, verbose:bool=False):
    output = path.join(output_loc, f'{db_name}_results.unprocessed.b6')
    run_process(['hmmsearch', '--domtblout', output, '--cpu', str(threads), db_loc, genes_faa], logger, verbose=verbose)
    # Parse hmmsearch output
    if not (path.isfile(output) and stat(output).st_size > 0):
        return pd.DataFrame()
    hits = parse_hmmsearch_domtblout(output)
    if len(hits) < 1:
        return pd.DataFrame()
    return formater(hits)


def get_sig_row(row, evalue_lim:float=1e-15):
    """Check if hmm match is significant, based on dbCAN described parameters"""
    tstart, tend, tlen, evalue = row[['target_start', 'target_end', 'target_length', 'full_evalue']].values
    perc_cov = (tend - tstart)/tlen
    if perc_cov >= .35 and evalue <= evalue_lim:
        return True
    else:
        return False


#TODO decide if we need use_hmmer_thresholds:bool=False
def generic_hmmscan_formater(hits:pd.DataFrame,  db_name:str, hmm_info_path:str=None, top_hit:bool=True):
    if hmm_info_path is None:
        hmm_info = None
        hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
        hits_sig = sig_scores(hits, hmm_info)
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    if top_hit:
        # Get the best hits
        hits_sig = hits_sig.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    hits_df = hits_sig[['target_id', 'query_id']]
    hits_df.set_index('query_id', inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns = [f"{db_name}_id"]
    if hmm_info is not None:
        hits_df = hits_df.merge(hmm_info[['definition']], how='left', left_on=f"{db_name}_id", right_index=True)
        hits_df.rename(columns={'definition': f"{db_name}_hits"}, inplace=True)
    return hits_df


def sig_scores(hits:pd.DataFrame, score_db:pd.DataFrame) -> pd.DataFrame:
    is_sig = list()
    for i, frame in hits.groupby('target_id'):
        row = score_db.loc[i]
        if row['score_type'] == 'domain':
            score = frame.domain_score
        elif row['score_type'] == 'full':
            score = frame.full_score
        elif row['score_type'] == '-':
            continue
        else:
            raise ValueError(row['score_type'])
        frame = frame.loc[score.astype(float) > float(row.threshold)]
        is_sig.append(frame)
    if len(is_sig) > 0:
        return pd.concat(is_sig)
    else:
        return pd.DataFrame()



