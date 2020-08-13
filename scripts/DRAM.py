#!/usr/bin/env python3

import argparse

from mag_annotator.annotate_bins import annotate_bins_cmd, annotate_called_genes
from mag_annotator.summarize_genomes import summarize_genomes
from mag_annotator.pull_sequences import pull_sequences, get_gene_neighborhoods

# TODO: refactor parses to limit duplication

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()
    annotate_parser = subparsers.add_parser('annotate', help="Annotate genomes/contigs/bins/MAGs",
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    annotate_genes_parser = subparsers.add_parser('annotate_genes', help="Annotate already called genes, limited "
                                                                         "functionality compared to annotate",
                                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    distill_parser = subparsers.add_parser('distill', help="Summarize metabolic content of annotated genomes",
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    strainer_parser = subparsers.add_parser('strainer', help="Strain annotations down to genes of interest",
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    neighborhood_parser = subparsers.add_parser('neighborhoods', help="Find neighborhoods around genes of interest",
                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser for annotating mags, you know the real thing
    annotate_parser.add_argument('-i', '--input_fasta',
                                 help="fasta file, optionally with wildcards to point to multiple fastas",
                                 required=True)
    annotate_parser.add_argument('-o', '--output_dir', help="output directory")
    annotate_parser.add_argument('--min_contig_size', type=int, default=2500,
                                 help='minimum contig size to be used for gene prediction')
    prodigal_mode_choices = ['train', 'meta', 'single']
    annotate_parser.add_argument('--prodigal_mode', type=str, default='meta', choices=prodigal_mode_choices,
                                 help='Mode of prodigal to use for gene calling. NOTE: normal or single mode require '
                                      'genomes which are high quality with low contamination and long contigs (average '
                                      'length >3 Kbp).')
    # prodigal_trans_table_choices = ['auto'] + [str(i) for i in range(1, 26)]
    prodigal_trans_table_choices = [str(i) for i in range(1, 26)]
    annotate_parser.add_argument('--trans_table', type=str, default='11', choices=prodigal_trans_table_choices,
                                 help='Translation table for prodigal to use for gene calling.')
    annotate_parser.add_argument('--bit_score_threshold', type=int, default=60,
                                 help='minimum bitScore of search to retain hits')
    annotate_parser.add_argument('--rbh_bit_score_threshold', type=int, default=350,
                                 help='minimum bitScore of reverse best hits to retain hits')
    annotate_parser.add_argument('--custom_db_name', action='append', help="Names of custom databases, can be used"
                                                                           "multiple times.")
    annotate_parser.add_argument('--custom_fasta_loc', action='append',
                                 help="Location of fastas to annotated against, can be used multiple times but"
                                      "must match nubmer of custom_db_name's")
    annotate_parser.add_argument('--gtdb_taxonomy', action='append', default=[],
                                 help='Summary file from gtdbtk taxonomy assignment from bins, can be used multiple'
                                      'times')
    annotate_parser.add_argument('--checkm_quality', action='append', default=[],
                                 help='Summary of of checkM quality assessment from bins, can be used multiple times')
    annotate_parser.add_argument('--use_uniref', action='store_true', default=False,
                                 help='Annotate these fastas against UniRef, drastically decreases run time and memory '
                                      'requirements')
    annotate_parser.add_argument('--low_mem_mode', action='store_true', default=False,
                                 help='Skip annotating with uniref and use kofam instead of KEGG genes even if '
                                      'provided. Drastically decreases memory usage')
    annotate_parser.add_argument('--skip_trnascan', action='store_true', default=False)
    annotate_parser.add_argument('--keep_tmp_dir', action='store_true', default=False)
    annotate_parser.add_argument('--threads', type=int, default=10, help='number of processors to use')
    annotate_parser.add_argument('--verbose', action='store_true', default=False)
    annotate_parser.set_defaults(func=annotate_bins_cmd)

    # parser for annotating already called genes
    annotate_genes_parser.add_argument('-i', '--input_faa', help="fasta file, optionally with wildcards to point to"
                                                                 "individual MAGs", required=True)
    annotate_genes_parser.add_argument('-o', '--output_dir', help="output directory")
    annotate_genes_parser.add_argument('--bit_score_threshold', type=int, default=60,
                                       help='minimum bitScore of search to retain hits')
    annotate_genes_parser.add_argument('--rbh_bit_score_threshold', type=int, default=350,
                                       help='minimum bitScore of reverse best hits to retain hits')
    annotate_genes_parser.add_argument('--custom_db_name', action='append', help="Names of custom databases, can be "
                                                                                 "used multiple times.")
    annotate_genes_parser.add_argument('--custom_fasta_loc', action='append',
                                       help="Location of fastas to annotated against, can be used multiple times but"
                                            "must match nubmer of custom_db_name's")
    annotate_genes_parser.add_argument('--use_uniref', action='store_true', default=False,
                                       help='Annotate these fastas against UniRef, drastically decreases run time and '
                                            'memory requirements')
    annotate_genes_parser.add_argument('--low_mem_mode', action='store_true', default=False,
                                       help='Skip annotating with uniref and use kofam instead of KEGG genes even if '
                                            'provided. Drastically decreases memory usage')
    annotate_genes_parser.add_argument('--keep_tmp_dir', action='store_true', default=False)
    annotate_genes_parser.add_argument('--threads', type=int, default=10, help='number of processors to use')
    annotate_genes_parser.add_argument('--verbose', action='store_true', default=False)
    annotate_genes_parser.set_defaults(func=annotate_called_genes)

    # parser for summarizing genomes
    distill_parser.add_argument("-i", "--input_file", help="Annotations path")
    distill_parser.add_argument("-o", "--output_dir", help="Directory to write summarized genomes")
    distill_parser.add_argument("--rrna_path", help="rRNA output from annotation")
    distill_parser.add_argument("--trna_path", help="tRNA output from annotation")
    distill_parser.add_argument("--groupby_column", help="Column from annotations to group as organism units",
                                default='fasta')
    distill_parser.set_defaults(func=summarize_genomes)

    # parser for getting genes
    input_group = strainer_parser.add_argument_group('Input and output files')
    input_group.add_argument('-i', '--input_annotations', required=True, help='annotations file to pull genes from')
    input_group.add_argument('-f', '--input_fasta', required=True, help='fasta file to filter')
    input_group.add_argument('-o', '--output_fasta', default='pull_genes.fasta',
                             help='location to write filtered fasta')
    name_group = strainer_parser.add_argument_group('Specific names to keep')
    name_group.add_argument('--fastas', nargs='*', default=None, help='space separated list of fastas to keep')
    name_group.add_argument('--scaffolds', nargs='*', default=None,
                            help='space separated list of scaffolds to keep')
    name_group.add_argument('--genes', nargs='*', default=None, help='space separated list of genes to keep')
    annotation_group = strainer_parser.add_argument_group('Annotation filters')
    annotation_group.add_argument('--identifiers', nargs='*', default=None, help='database identifiers to keep')
    annotation_group.add_argument('--categories', nargs='*', default=None,
                                  help='distillate categories to keep genes from')
    dram_group = strainer_parser.add_argument_group('DRAM based filters')
    dram_group.add_argument('--taxonomy', nargs='*', default=None,
                            help='Level of GTDBTk taxonomy to keep (e.g. c__Clostridia), space separated list')
    dram_group.add_argument('--completeness', default=None, type=float,
                            help='Minimum completeness of genome to keep genes')
    dram_group.add_argument('--contamination', default=None, type=float,
                            help='Maximum contamination of genome to keep genes')
    strainer_parser.set_defaults(func=pull_sequences)

    # parser for getting gene neighborhoods
    neighborhood_parser.add_argument("-i", "--input_file", help="Annotations path")
    neighborhood_parser.add_argument("-o", "--output_dir", help="Directory to write gene neighborhoods")
    neighborhood_parser.add_argument("--genes", nargs='*', help="Gene names from DRAM to find neighborhoods around")
    neighborhood_parser.add_argument("--identifiers", nargs='*',
                                     help="Database identifiers assigned by DRAM to find neighborhoods around")
    neighborhood_parser.add_argument("--categories", help="Distillate categories to build gene neighborhoods around.")
    neighborhood_parser.add_argument("--genes_loc", help="Location of genes.fna/genes.faa file to filter to "
                                                         "neighborhoods")
    neighborhood_parser.add_argument("--scaffolds_loc", help="Location of scaffolds.fna file to filter to "
                                                             "neighborhoods")
    neighborhood_parser.add_argument("--distance_genes", type=int, help="Number of genes away from center to include "
                                                                        "in neighborhoods")
    neighborhood_parser.add_argument("--distance_bp", type=int, help="Number of genes away from center to include "
                                                                     "in neighborhoods")
    neighborhood_parser.set_defaults(func=get_gene_neighborhoods)

    args = parser.parse_args()
    args_dict = {i: j for i, j in vars(args).items() if i != 'func'}
    args.func(**args_dict)
