import re
import subprocess
from collections import Counter
from glob import glob
from os import path
from pkg_resources import resource_filename
import json
from urllib.request import urlopen
import warnings


def get_config_loc():
    return path.abspath(resource_filename('mag_annotator', 'CONFIG'))


def get_database_locs(config_loc=None):
    if config_loc is None:
        config_loc = get_config_loc()
    return json.loads(open(config_loc).read())


def download_file(url, output_file=None, verbose=True):
    if verbose:
        print('downloading %s' % url)
    if output_file is None:
        return urlopen(url).read().decode('utf-8')
    else:
        run_process(['wget', '-O', output_file, url], verbose=verbose)


def run_process(command, shell=False, capture_stdout=True, check=True, verbose=False):
    """Standardization of parameters for using subprocess.run, provides verbose mode and option to run via shell"""
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    if capture_stdout:
        return subprocess.run(command, check=check, shell=shell, stdout=subprocess.PIPE,
                              stderr=stderr).stdout.decode(errors='ignore')
    else:
        subprocess.run(command, check=check, shell=shell, stdout=stdout, stderr=stderr)


def make_mmseqs_db(fasta_loc, output_loc, create_index=True, threads=10, verbose=False):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2"""
    run_process(['mmseqs', 'createdb', fasta_loc, output_loc], verbose=verbose)
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), 'tmp')
        run_process(['mmseqs', 'createindex', output_loc, tmp_dir, '--threads', str(threads)], verbose=verbose)


def multigrep(search_terms, search_against, split_char='\n', output='.'):
    # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, 'hits.txt')
    with open(hits_file, 'w') as f:
        f.write('%s\n' % '\n'.join(search_terms))
    results = run_process(['grep', '-a', '-F', '-f', hits_file, search_against], capture_stdout=True, verbose=False)
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


def get_ids_from_annotation(frame):
    id_list = list()
    # get kegg ids
    if 'kegg_id' in frame:
        id_list += [j.strip() for i in frame.kegg_id.dropna() for j in i.split(',')]
    # get kegg ec numbers
    if 'kegg_hit' in frame:
        for kegg_hit in frame.kegg_hit.dropna():
            id_list += [i[1:-1] for i in re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', kegg_hit)]
    # get merops ids
    if 'peptidase_family' in frame:
        id_list += [j.strip() for i in frame.peptidase_family.dropna() for j in i.split(';')]
    # get cazy ids
    if 'cazy_hits' in frame:
        id_list += [j.split(' ')[0] for i in frame.cazy_hits.dropna() for j in i.split(';')]
        # get cazy ec numbers
        for cazy_hit in frame.cazy_hits.dropna():
            id_list += [i[1:-1].replace(' ', ':') for i in re.findall(r'\(EC \d*.\d*.\d*.\d*\)', cazy_hit)]
    # get pfam ids
    if 'pfam_hits' in frame:
        id_list += [j[1:-1].split('.')[0] for i in frame.pfam_hits.dropna()
                    for j in re.findall(r'\[PF\d\d\d\d\d.\d*\]', i)]
    return Counter(id_list)


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
        return text[:len(suffix)]
    return text  # or whatever
