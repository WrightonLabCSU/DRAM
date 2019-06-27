import argparse

from mag_annotator.summarize_genomes import summarize_genomes

genome_summary_frame_path = '/Users/shafferm/lab/AMG/genome_summary_table.tsv'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input_file", help="Annotations path")
    parser.add_argument("-o", "--output_dir", help="Directory to write summarized genomes")
    parser.add_argument("--genome_summary_path", help="Frame of genome summary to be filled in",
                        default=genome_summary_frame_path)
    parser.add_argument("--trna_path", help="tRNA output from annotation")
    parser.add_argument("--modules_path", help="File outlining metabolisms of modules from KEGG")
    parser.add_argument("--group_column", help="Column from annotations to group as organism units", default='fasta')
    parser.add_argument("--viral", default=False, action='store_true',
                        help="If sample is viral will remove empty functions")
    parser.add_argument("--min_cov", type=float, default=.001, help="Minimum coverage to include module in summary")

    args = parser.parse_args()
    summarize_genomes(args.input_file, args.genome_summary_frame, args.trna_path, args.modules_path, args.group_column,
                      args.output_file, args.viral, args.min_cov)
