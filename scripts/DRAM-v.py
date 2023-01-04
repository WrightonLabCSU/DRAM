#!/usr/bin/env python3

import argparse

from mag_annotator.annotate_vgfs import annotate_vgfs, remove_bad_chars
from mag_annotator.summarize_vgfs import summarize_vgfs
from mag_annotator.pull_sequences import pull_sequences, get_gene_neighborhoods

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()
    annotate_parser = subparsers.add_parser('annotate', help="Annotate viral contigs",
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    distill_parser = subparsers.add_parser('distill',
                                           help="Summarize AMGs in annotated viral genome fragments",
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    strainer_parser = subparsers.add_parser('strainer', help="Strain annotations down to genes of interest",
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    remove_parser = subparsers.add_parser('remove_bad_characters', help="Removes ; and = from fasta headers and "
                                                                        "VIRSorter_affi-contigs.tab files",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    neighborhood_parser = subparsers.add_parser('neighborhoods', help="Find neighborhoods around genes of interest",
                                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser for annotating vgfs
    annotate_parser.add_argument('-i', '--input_fasta', help="fasta file, output from ", required=True)
    annotate_parser.add_argument('-v', '--virsorter_affi_contigs', help="VirSorter VIRSorter_affi-contigs.tab "
                                                                        "output file")
    annotate_parser.add_argument('-o', '--output_dir', help="output directory")
    annotate_parser.add_argument('--min_contig_size', type=int, default=2500,
                                 help='minimum contig size to be used for gene prediction')
    annotate_parser.add_argument('--split_contigs', action='store_true', default=False,
                                 help='Split contigs from input fasta into separate')
    prodigal_mode_choices = ['train', 'meta', 'single']
    annotate_parser.add_argument('--prodigal_mode', type=str, default='meta', choices=prodigal_mode_choices,
                                 help='Mode of prodigal to use for gene calling. NOTE: normal or single mode require '
                                      'genomes which are high quality with low contamination and long contigs (average '
                                      'length >3 Kbp).')
    # prodigal_trans_table_choices = ['auto'] + [str(i) for i in range(1, 26)]
    prodigal_trans_table_choices = [str(i) for i in range(1, 26)]
    annotate_parser.add_argument('--trans_table', type=str, default='11', choices=prodigal_trans_table_choices,
                                 help='Translation table for prodigal to use for gene calling')
    annotate_parser.add_argument('--bit_score_threshold', type=int, default=60,
                                 help='minimum bitScore of search to retain hits')
    annotate_parser.add_argument('--rbh_bit_score_threshold', type=int, default=350,
                                 help='minimum bitScore of reverse best hits to retain hits')
    annotate_parser.add_argument('--kofam_use_dbcan2_thresholds', action='store_true', default=False,
                                 help='Use dbcan2 suggested HMM cutoffs for KOfam annotation instead of KOfam '
                                      'recommended cutoffs. This will be ignored if annotating with KEGG Genes.')
    annotate_parser.add_argument('--custom_db_name', action='append', default=[], # empty list is not bug, can't be changed
                                 help="Names of custom databases, can be used multiple times.")
    annotate_parser.add_argument('--custom_fasta_loc', action='append', default=[], # empty list is not bug, can't be changed
                                 help="Location of fastas to annotate against, can be used multiple times but"
                                      "must match nubmer of custom_db_name's")
    annotate_parser.add_argument('--custom_hmm_name', action='append',  default=[], # empty list is not bug, can't be changed
                                 help="Names of custom hmm databases, can be used multiple times.")
    annotate_parser.add_argument('--custom_hmm_loc', action='append',  default=[], # empty list is not bug, can't be changed
                                 help="Location of hmms to annotate against, can be used multiple times but"
                                      "must match nubmer of custom_hmm_name's")
    annotate_parser.add_argument('--custom_hmm_cutoffs_loc', action='append',
                                       help="Location of file with custom HMM cutoffs and descriptions, can be used "
                                            "multiple times.")
    annotate_parser.add_argument('--use_uniref', action='store_true', default=False,
                                 help='Annotate these fastas against UniRef, drastically increases'
                                      ' run time and memory requirements')
    annotate_parser.add_argument('--low_mem_mode', action='store_true', default=False,
                                 help='Skip annotating with uniref and use kofam instead of KEGG genes even if '
                                      'provided. Drastically decreases memory usage')
    annotate_parser.add_argument('--skip_trnascan', action='store_true', default=False)
    annotate_parser.add_argument('--keep_tmp_dir', action='store_true', default=False)
    annotate_parser.add_argument('--threads', type=int, default=10, help='number of processors to use')
    annotate_parser.add_argument('--verbose', action='store_true', default=False)
    annotate_parser.add_argument('--config_loc', default=None,
                                 help='location of an alternive config file that will over write the original at run time,'
                                 ' but not be saved or modified')
    annotate_parser.set_defaults(func=annotate_vgfs)

    # parser for summarizing genomes
    distill_parser.add_argument("-i", "--input_file", help="Annotations path")
    distill_parser.add_argument("-o", "--output_dir", help="Directory to write summarized genomes")
    distill_parser.add_argument("--groupby_column", help="Column from annotations to group as VGF units",
                                default='scaffold')
    distill_parser.add_argument("--max_auxiliary_score", type=int, default=3,
                                help="Maximum auxiliary score to consider gene as potential AMG")
    distill_parser.add_argument("--remove_transposons", default=False, action='store_true',
                                help="Do not consider genes on scaffolds with transposons as potential AMGs")
    distill_parser.add_argument("--remove_fs", default=False, action='store_true',
                                help="Do not consider genes near ends of scaffolds as potential AMGs")
    distill_parser.add_argument('--log_file_path', 
                                       help="A name and loctation for the log file")
    # distill_parser.add_argument("--remove_js", default=False, action='store_true',
    #                             help="Do not consider genes on possible non-viral contigs as potential AMGs")
    distill_parser.add_argument("--custom_distillate", help="Custom distillate form to add your own modules")
    distill_parser.add_argument('--config_loc', default=None,
                                 help='location of an alternive config file that will over write the original at run time,'
                                 ' but not be saved or modified')
    distill_parser.set_defaults(func=summarize_vgfs)

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
    dramv_group = strainer_parser.add_argument_group('DRAM-v based filters')
    dramv_group.add_argument('--amg_flags', default=None,
                             help='AMG flags to keep, if any one is present then it will be kept')
    dramv_group.add_argument('--aux_scores', nargs='*', default=None, type=int,
                             help='space separate list of auxiliary scores to keep')
    dramv_group.add_argument('--virsorter_category', nargs='*', default=None,
                             help='space separate list of virsorter categories to keep')
    amg_group = strainer_parser.add_argument_group('AMG filtering')
    amg_group.add_argument('-a', '--putative_amgs', default=False, action='store_true',
                           help='Only keep genes considered putative AMGs')
    amg_group.add_argument('--max_auxiliary_score', type=int, default=3,
                           help="Maximum auxiliary score to consider gene as potential AMG")
    amg_group.add_argument('--remove_transposons', default=False, action='store_true',
                           help="Do not consider genes on scaffolds with transposons as potential AMGs")
    amg_group.add_argument('--remove_fs', default=False, action='store_true',
                           help="Do not consider genes near ends of scaffolds as potential AMGs")
    # amg_group.add_argument("--remove_js", default=False, action='store_true',
    #                        help="Do not consider genes on possible non-viral contigs as potential AMGs")
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

    fasta_or_affi = remove_parser.add_mutually_exclusive_group(required=True)
    fasta_or_affi.add_argument('-i', '--input_fasta', help='Fasta file to remove ; and = from headers')
    fasta_or_affi.add_argument('-v', '--input_virsorter_affi_contigs', help='Fasta file to remove ; and = from headers')
    remove_parser.add_argument('-o', '--output', required=True, help='Name of output file. If cleaning a fasta file the'
                                                                     ' output file name must have no = or ;.')
    remove_parser.set_defaults(func=remove_bad_chars)

    args = parser.parse_args()
    args_dict = {i: j for i, j in vars(args).items() if i != 'func'}
    args.func(**args_dict)
