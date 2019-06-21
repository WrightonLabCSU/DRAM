import subprocess
from glob import glob
from os import path, remove


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


def multigrep(search_terms, search_against, output='.'):  # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, 'hits.txt')
    with open(hits_file, 'w') as f:
        f.write('%s\n' % '\n'.join(search_terms))
    results = subprocess.run(['grep', '-a', '-F', '-f', hits_file, search_against], stdout=subprocess.PIPE, check=True)
    processed_results = [i.strip()[1:] for i in results.stdout.decode('ascii').split('\n')]
    remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ''}


def merge_files(files_to_merge, outfile):
    """It's in the name"""
    with open(outfile, 'w') as outfile_handle:
        for file in glob(files_to_merge):
            with open(file) as f:
                outfile_handle.write(f.read())


def merge_files_w_header(gtf_files, outfile):
    """Merge files but keep the header from the first file, assumes files all have same header"""
    gtf_files = glob(gtf_files)
    with open(outfile, 'w') as f:
        f.write(open(gtf_files[0]).readline())
        for gtf in gtf_files:
            content = ''.join(open(gtf).readlines()[1:])
            f.write(content)
