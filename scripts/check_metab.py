import argparse

from checkMetab.check_metab import main


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--annotation_file', help="annotation file")
    parser.add_argument('-m', '--metab_database', help="metabolism form")
    parser.add_argument('-o', '--output_file', default='modules_present.tsv', help='Output file')
    parser.add_argument('-c', '--min_cov', default=0.5, type=int, help='Minimum coverage value for o')

    args = parser.parse_args()

    main(args.annotation_file, args.metab_database, args.output_file, args.min_cov)
