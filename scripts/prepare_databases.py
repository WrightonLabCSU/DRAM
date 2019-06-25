import argparse

from checkMetab.database_processing import prepare_databases


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output_dir', default="~/checkMetab_data", help="output directory")
    parser.add_argument('--kegg_loc', default=None,
                        help="KEGG protein file, should be a single .pep, please merge all KEGG pep files")
    parser.add_argument('--kegg_download_date', default=None, help="Date KEGG was download to include in database name")
    parser.add_argument('--uniref_loc', default=None, help="File path to uniref, if already downloaded")
    parser.add_argument('--uniref_version', default='90', help="UniRef version to download")
    parser.add_argument('--pfam_loc', default=None, help="File path to pfam-A hmm file, if already downloaded")
    parser.add_argument('--pfam_release', default='32.0', help="Pfam release to download")
    parser.add_argument('--dbcan_loc', default=None, help="File path to dbCAN, if already downloaded")
    parser.add_argument('--dbcan_version', default='7', type=str, help='version of dbCAN to use')
    parser.add_argument('--viral_loc', default=None, help="File path to viral protein faa, if already downloaded")
    parser.add_argument('--keep_db_files', default=False, action='store_true', help="Keep unporcessed database files")
    parser.add_argument('--threads', default=10, type=int, help="Number of threads to use building mmseqs2 databases")
    parser.add_argument('--verbose', default=False, action='store_true', help="Make it talk more")

    args = parser.parse_args()

    output_dir = args.output_dir
    kegg_loc = args.kegg_loc
    kegg_download_date = args.kegg_download_date
    uniref_loc = args.uniref_loc
    uniref_version = args.uniref_version
    pfam_loc = args.pfam_loc
    pfam_release = args.pfam_release
    dbcan_loc = args.dbcan_loc
    dbcan_version = args.dbcan_version
    viral_loc = args.viral_loc
    keep_db_files = args.keep_db_files
    threads = args.threads
    verbose = args.verbose

    prepare_databases(output_dir=output_dir, kegg_loc=kegg_loc, kegg_download_date=kegg_download_date,
                      uniref_loc=uniref_loc, uniref_version=uniref_version, pfam_loc=pfam_loc,
                      pfam_version=pfam_release, dbcan_loc=dbcan_loc, dbcan_version=dbcan_version, viral_loc=viral_loc,
                      keep_database_files=keep_db_files, threads=threads, verbose=verbose)
