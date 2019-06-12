from os import path
from checkMetab.annotate_scaffolds import run_process


def download_unifref(output_dir, uniref_version='90', verbose=True):
    """"""
    if verbose:
        print('downloading uniref fasta to %s' % output_dir)
    uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % uniref_version)
    run_process(['wget', '-O', uniref_fasta_zipped,
                 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz'
                 % (uniref_version, uniref_version)], verbose=verbose)
    if verbose:
        print('unzipping %s' % uniref_fasta_zipped)
    run_process(['gunzip', uniref_fasta_zipped], verbose=verbose)


def download_and_process_pfam(output_dir, pfam_release='32.0', threads=10, verbose=True):
    if verbose:
        print('downloading pfam msa to %s' % output_dir)
    pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
    run_process(['wget', '-O', pfam_full_zipped,
                 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam%s/Pfam-A.full.gz' % pfam_release],
                verbose=verbose)
    mmseq_msa = path.join(output_dir, 'pfam.mmsmsa')
    run_process(['mmseqs', 'convertmsa', pfam_full_zipped, mmseq_msa], verbose=verbose)
    mmseq_profile = path.join(output_dir, 'pfam.mmspro')
    run_process(['mmseqs', 'msa2profile', mmseq_msa, mmseq_profile, '--match-mode', '1', '--threads', str(threads)],
                verbose=verbose)
    tmp_dir = path.join(output_dir, 'tmp')
    run_process(['mmseqs', 'createindex', mmseq_profile, tmp_dir, '-k', '5', '-s', '7', '--threads', str(threads)],
                verbose=verbose)
