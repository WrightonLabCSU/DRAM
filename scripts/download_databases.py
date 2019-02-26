import argparse

from checkMetab.annotate_scaffolds import download_and_process_pfam, download_unifref


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_dir', default='.', help="output directory")
    parser.add_argument('--pfam_release', default='32.0', help='Pfam release to download')
    parser.add_argument('--uniref_version', default='90', help='UniRef version to download')

    args = parser.parse_args()

    fasta_loc = args.input_fasta
    output = args.output_dir
    min_contig_size = args.min_contig_size

    download_and_process_pfam(args.output_dir, pfam_release=args.pfam_release)
    download_unifref(args.output_dir, uniref_version=args.uniref_version)
