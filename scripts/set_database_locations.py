import argparse

from mag_annotator.database_processing import set_database_paths

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--kegg_loc', default=None, help='mmseqs2 database file from kegg .pep file')
    parser.add_argument('--uniref_loc', default=None, help='mmseqs2 database file from uniref .faa')
    parser.add_argument('--pfam_loc', default=None, help='mmseqs2 database file from pfam .hmm')
    parser.add_argument('--dbcan_loc', default=None, help='hmm file for dbcan, already processed with hmmpress')
    parser.add_argument('--viral_loc', default=None, help='mmseqs2 database file from ref seq viral gene collection')

    args = parser.parse_args()

    kegg_loc = args.kegg_loc
    uniref_loc = args.uniref_loc
    pfam_loc = args.pfam_loc
    dbcan_loc = args.dbcan_loc
    viral_loc = args.viral_loc

    set_database_paths(kegg_loc, uniref_loc, pfam_loc, dbcan_loc, viral_loc)
