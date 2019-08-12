import subprocess
from glob import glob
from os import path, remove
from pkg_resources import resource_filename
import json
from urllib.request import urlopen
import warnings


def get_database_locs():
    return json.loads(open(path.abspath(resource_filename('mag_annotator', 'CONFIG'))).read())


def download_file(url, output_file=None, verbose=True):
    if verbose:
        print('downloading %s' % url)
    if output_file is None:
        return urlopen(url).read().decode('utf-8')
    else:
        run_process(['wget', '-O', output_file, url], verbose=verbose)


def run_process(command, shell=False, verbose=False):
    """Standardization of parameters for using subprocess.run, provides verbose mode and option to run via shell"""
    stdout = subprocess.DEVNULL
    stderr = subprocess.DEVNULL
    if verbose:
        stdout = None
        stderr = None
    subprocess.run(command, check=True, shell=shell, stdout=stdout, stderr=stderr)


def make_mmseqs_db(fasta_loc, output_loc, create_index=True, threads=10, verbose=False):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2"""
    run_process(['mmseqs', 'createdb', fasta_loc, output_loc], verbose=verbose)
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), 'tmp')
        run_process(['mmseqs', 'createindex', output_loc, tmp_dir, '--threads', str(threads)], verbose=verbose)


def multigrep(search_terms, search_against, output='.', skip_chars=1):
    # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, 'hits.txt')
    with open(hits_file, 'w') as f:
        f.write('%s\n' % '\n'.join(search_terms))
    results = subprocess.run(['grep', '-a', '-F', '-f', hits_file, search_against], stdout=subprocess.PIPE, check=True)
    processed_results = [header.strip()[skip_chars:] if i > 0 else header.strip() for i, header in
                         enumerate(results.stdout.decode('ascii').split('\n'))]
    processed_results = [i.strip() for i in results.stdout.decode(errors='ignore').strip().split(' \n\x00')
                         if len(i) > 0]
    # remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ''}


def merge_files(paths_to_files_to_merge, outfile, has_header=False):
    """It's in the name, if has_header assumes all files have the same header"""
    files_to_merge = glob(paths_to_files_to_merge)
    if len(files_to_merge) == 0:
        warnings.warn('No files to merge found with path %s' % paths_to_files_to_merge)
    else:
        with open(outfile, 'w') as outfile_handle:
            if has_header:
                outfile_handle.write(open(files_to_merge[0]).readline())
            for file in files_to_merge:
                with open(file) as f:
                    if has_header:
                        _ = f.readline()
                    outfile_handle.write(f.read())
