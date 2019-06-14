from os import path
from checkMetab.annotate_bins import run_process


def download_file(url, output_file, verbose=True):
    if verbose:
        print('downloading %s' % url)
    run_process(['wget', '-O', output_file, url], verbose=verbose)


def download_unifref(output_dir, uniref_version='90', verbose=True):
    """"""
    uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
    uniref_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' %\
                 (uniref_version, uniref_version)
    download_file(uniref_url, uniref_fasta_zipped)
    if verbose:
        print('unzipping %s' % uniref_fasta_zipped)
    run_process(['gunzip', uniref_fasta_zipped], verbose=verbose)


def process_mmspro(full_alignment, output_dir, db_name='db', threads=10, verbose=True):
    mmseq_msa = path.join(output_dir, '%s.mmsmsa' % db_name)
    run_process(['mmseqs', 'convertmsa', full_alignment, mmseq_msa], verbose=verbose)
    mmseq_profile = path.join(output_dir, '%s.mmspro' % db_name)
    run_process(['mmseqs', 'msa2profile', mmseq_msa, mmseq_profile, '--match-mode', '1', '--threads', str(threads)],
                verbose=verbose)
    tmp_dir = path.join(output_dir, 'tmp')
    run_process(['mmseqs', 'createindex', mmseq_profile, tmp_dir, '-k', '5', '-s', '7', '--threads', str(threads)],
                verbose=verbose)
    return mmseq_profile


def download_and_process_pfam(output_dir, pfam_release='32.0', threads=10, verbose=True):
    pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
    download_file('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.full.gz' % pfam_release, pfam_full_zipped)
    process_mmspro(pfam_full_zipped, output_dir, 'pfam', threads, verbose)
