import re
import subprocess
from collections import Counter
from os import path
from pkg_resources import resource_filename
import json
from urllib.request import urlopen
import pandas as pd


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


# unify this with get_ids_from_annotation
def get_ids_from_row(row):
    id_list = list()
    # get kegg ids
    if 'kegg_id' in row and not pd.isna(row['kegg_id']):
        id_list += [j for j in row['kegg_id'].split(',')]
    # get ec numbers
    if 'kegg_hit' in row and not pd.isna(row['kegg_hit']):
        id_list += [i[1:-1] for i in re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', row['kegg_hit'])]
    # get merops ids
    if 'peptidase_family' in row and not pd.isna(row['peptidase_family']):
        id_list += [j for j in row['peptidase_family'].split(';')]
    # get cazy ids
    if 'cazy_hits' in row and not pd.isna(row['cazy_hits']):
        id_list += [j.strip().split(' ')[0] for j in row['cazy_hits'].split(';')]
    return set(id_list)


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


def get_genes_from_identifiers(annotations, genes=None, fastas=None, scaffolds=None, identifiers=None, categories=None):
    specific_genes_to_keep = list()
    # filter fastas
    if fastas is not None:
        for fasta in fastas:
            specific_genes_to_keep += list(annotations.loc[annotations['fasta'] == fasta].index)
    # filter scaffolds
    # TODO: remove this functionality or modify since scaffolds are guaranteed unique
    if scaffolds is not None:
        for scaffold in scaffolds:
            specific_genes_to_keep += list(annotations.loc[annotations['scaffold'] == scaffold].index)
    # filter genes
    if genes is not None:
        specific_genes_to_keep += genes
    # filter down annotations based on specific genes
    if len(specific_genes_to_keep) > 0:
        annotations_to_keep = annotations.loc[specific_genes_to_keep]
    else:
        annotations_to_keep = annotations

    # filter based on annotations
    if (identifiers is not None) or (categories is not None):
        annotation_genes_to_keep = list()

        # make a dictionary of genes to
        gene_to_ids = dict()
        for i, row in annotations_to_keep.iterrows():
            row_ids = get_ids_from_row(row)
            if len(row_ids) > 0:
                gene_to_ids[i] = set(row_ids)

        # get genes with ids
        if identifiers is not None:
            identifiers = set(identifiers)
            for gene, ids in gene_to_ids.items():
                if len(set(ids) & set(identifiers)) > 0:
                    annotation_genes_to_keep.append(gene)

        # get genes from distillate categories
        if categories is not None:
            db_locs = get_database_locs()
            genome_summary_form = pd.read_csv(db_locs['genome_summary_form'], sep='\t')
            for level in ['module', 'sheet', 'header', 'subheader']:
                for category, frame in genome_summary_form.loc[~pd.isna(genome_summary_form[level])].groupby(level):
                    if category in categories:
                        for gene, ids in gene_to_ids.items():
                            if len(ids & set(frame['gene_id'])) > 0:
                                annotation_genes_to_keep.append(gene)
    else:
        annotation_genes_to_keep = list(annotations_to_keep.index)
    return annotation_genes_to_keep
