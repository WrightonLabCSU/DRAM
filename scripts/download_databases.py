import argparse

from checkMetab.database_processing import prepare_databases


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_dir', default='.', help="output directory")
    parser.add_argument('--kegg_loc', default=None,
                        help="KEGG protein file, should be a single .pep, please merge all KEGG pep files")
    parser.add_argument('--kegg_download_date', default=None, help="Date KEGG was download to include in database name")
    parser.add_argument('--pfam_release', default='32.0', help="Pfam release to download")
    parser.add_argument('--uniref_version', default='90', help="UniRef version to download")
    parser.add_argument('--keep_db_files', default=False, action='store_true', help="Keep unporcessed database files")
    parser.add_argument('--threads', default=10, type=int, help="Number of threads to use building mmseqs2 databases")
    parser.add_argument('--verbose', default=False, action='store_true', help="Make it talk more")

    args = parser.parse_args()

    output_dir = args.output_dir
    kegg_loc = args.kegg_loc
    kegg_download_date = args.kegg_download_date
    keep_db_files = args.keep_db_files
    threads = args.threads
    verbose = args.verbose

    prepare_databases(output_dir, kegg_loc, kegg_download_date, keep_db_files, threads, verbose)
