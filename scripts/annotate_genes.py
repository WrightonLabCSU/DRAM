#!/usr/bin/env python3

import argparse

from checkMetab.annotate_scaffolds import main


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_fasta', help="fasta file with scaffolds, should be nucleotides", required=True)
    parser.add_argument('--kegg_loc', help='mmseqs2 database file from kegg .pep file')
    parser.add_argument('--uniref_loc', help='mmseqs2 database file from uniref .faa')
    parser.add_argument('--pfam_loc', help='mmseqs2 database file from pfam .hmm')
    parser.add_argument('--dbcan_loc', help='hmm file for dbcan, already processed with hmmpress')
    parser.add_argument('--viral_loc', help='mmseqs2 database file from ref seq viral gene collection')
    parser.add_argument('-o', '--output_dir', help="output directory")
    parser.add_argument('--min_contig_size', type=int, default=5000,
                        help='minimum contig size to be used for gene prediction')
    parser.add_argument('--min_bitscore', type=int, default=60, help='minimum bitScore of search to retain hits')
    parser.add_argument('--strict_bitscore', type=int, default=350,
                        help='minimum bitScore of reverse best hits to retain hits')
    parser.add_argument('--keep_working_dir', action='store_true', default=False)
    parser.add_argument('--threads', type=int, default=10, help='number of processors to use')
    parser.add_argument('--verbose', action='store_true', default=False)

    args = parser.parse_args()

    fasta_loc = args.input_fasta
    kegg_loc = args.kegg_loc
    uniref_loc = args.uniref_loc
    pfam_loc = args.pfam_loc
    dbcan_loc = args.dbcan_loc
    viral_loc = args.viral_loc
    output = args.output_dir
    min_contig_size = args.min_contig_size
    min_bitscore = args.min_bitscore
    strict_bitscore = args.strict_bitscore
    keep_tmp = args.keep_working_dir
    threads = args.threads
    verbose = args.verbose

    main(fasta_loc, kegg_loc, uniref_loc, pfam_loc, dbcan_loc, viral_loc, output, min_contig_size, min_bitscore, strict_bitscore,
         keep_tmp, threads, verbose)
