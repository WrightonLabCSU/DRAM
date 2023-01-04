"""
Main control point for the viral annotation process
"""
from datetime import datetime
from os import path, mkdir
import re
import warnings
import logging
import pandas as pd
import numpy as np
from skbio.io import read as read_sequence
from skbio.io import write as write_sequence

from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.annotate_bins import annotate_fastas, get_fasta_name, perform_fasta_checks
from mag_annotator.utils import setup_logger
from mag_annotator.summarize_genomes import get_ids_from_annotations_all

VMAG_DBS_TO_ANNOTATE = ('kegg', 'kofam_hmm', 'kofam_ko_list', 'uniref', 'peptidase', 'pfam', 'dbcan', 'viral', 'vogdb')
VIRSORTER_COLUMN_NAMES = ['gene_name', 'start_position', 'end_position', 'length', 'strandedness',
                          'viral_protein_cluster_hit', 'viral_protein_cluster_hit_score',
                          'viral_protein_cluster_hit_evalue', 'viral_protein_cluster_category', 'pfam_hit',
                          'pfam_hit_score', 'pfam_hit_evalue', 'name']
VIRSORTER_HALLMARK_GENE_CATEGORIES = {'0', '3'}
VIRSORTER_VIRAL_LIKE_GENE_CATEGORIES = {'1', '4'}
TRANSPOSON_PFAMS = {'PF01609', 'PF00872', 'PF01610', 'PF01527', 'PF02371', 'PF01710', 'PF01385', 'PF01548', 'PF01526',
                    'PF01797', 'PF02899', 'PF05717', 'PF07592', 'PF03050', 'PF04754', 'PF04986', 'PF03400'}
CELL_ENTRY_CAZYS = {'CBM50', 'GH102', 'GH103', 'GH104', 'GH108', 'GH18', 'GH19', 'GH22', 'GH23', 'GH24', 'GH25', 'GH73',
                    'PL9', 'CBM12', 'CBM14', 'CBM18', 'CBM19', 'AA15', 'AA10', 'GH32', 'GH58', 'PL16'}
VIRAL_PEPTIDASES_MEROPS = {'A02H', 'A02G', 'A02F', 'A02E', 'A02D', 'A02C', 'A02B', 'A02A', 'A03B', 'A03A', 'A11B',
                           'A11A', 'A22B', 'A22A', 'A33', 'C01B', 'C01A', 'C02B', 'C02A', 'C03H', 'C03G', 'C03F',
                           'C03E', 'C03D', 'C03C', 'C03B', 'C03A', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C104',
                           'C105', 'C107', 'C108', 'C113', 'C14B', 'C14A', 'C16B', 'C16A', 'C18', 'C19', 'C21', 'C23',
                           'C24', 'C26', 'C27', 'C28', 'C30', 'C31', 'C32', 'C33', 'C36', 'C37', 'C39', 'C40', 'C42',
                           'C44', 'C46', 'C48', 'C49', 'C51', 'C53', 'C57', 'C59', 'C60B', 'C60A', 'C62', 'C63', 'C71',
                           'C74', 'C76', 'C85B', 'C85A', 'C87', 'C89', 'C97', 'C99', 'G02', 'I02', 'I04', 'I08', 'I24',
                           'I25C', 'I25B', 'I25A', 'I29', 'I32', 'I36', 'I43', 'I50B', 'I50A', 'I51', 'I63', 'I75',
                           'I87', 'I91', 'M03C', 'M03B', 'M03A', 'M10C', 'M10B', 'M10A', 'M12C', 'M12B', 'M12A', 'M13',
                           'M14C', 'M14D', 'M14B', 'M14A', 'M15D', 'M15C', 'M15B', 'M15A', 'M16C', 'M16B', 'M16A',
                           'M20E', 'M20A', 'M20D', 'M20F', 'M20B', 'M20C', 'M23B', 'M23A', 'M27', 'M38', 'M41', 'M42',
                           'M43B', 'M43A', 'M44', 'M48C', 'M48B', 'M48A', 'M56', 'M60', 'M67C', 'M67B', 'M67A', 'M78',
                           'M79', 'M86', 'N01', 'N02', 'N04', 'N05', 'N07', 'N08', 'N09', 'N10E', 'N10D', 'N10C',
                           'N10B', 'N11', 'S01F', 'S01E', 'S01D', 'S01C', 'S01B', 'S01A', 'S03', 'S06', 'S07', 'S08C',
                           'S08B', 'S08A', 'S09D', 'S09C', 'S09B', 'S09A', 'S11', 'S12', 'S14', 'S16', 'S21', 'S24',
                           'S26C', 'S26B', 'S26A', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33', 'S49C', 'S49B', 'S49A',
                           'S50', 'S53', 'S54', 'S62', 'S65', 'S69', 'S73', 'S74', 'S75', 'S77', 'S78', 'S80', 'S81',
                           'T01B', 'T01A', 'T03', 'U32', 'U40'}


def is_affi_tab_not_fasta(input_file: str) -> bool:
    """
    Checks that a file is a VIRSorter affi tab file.

    :param input_file:
        the input file that the user provided.
    :returns:
        True if there are 11 delimiters on the second line.
    """
    with open(input_file) as input_io:
        num_delimiter = input_io.readlines()[1].count('|')
    if num_delimiter == 11:
        return True
    return False


def remove_bad_chars_fasta(fasta):
    if is_affi_tab_not_fasta(fasta):
        raise ValueError("The input file format matches virsorter affi contigs not fasta, please check "
                         "that the flags match the file type: '-v' for virsorter, and '-i' for fasta.")

    new_headers = list()
    new_seqs = list()
    for seq in read_sequence(fasta, format='fasta'):
        new_header = seq.metadata['id']
        new_header = new_header.replace(';', '_').replace('=', '_')
        if new_header in new_headers:
            raise ValueError('Replacement of ; or = with _ generated redundant sequence headers, must be modified '
                             'manually')
        seq.metadata['id'] = new_header
        new_seqs.append(seq)
    return new_seqs


def remove_bad_chars_virsorter_affi_contigs(virsorter_in: str) -> str:
    """
    Determine and call the appropriate function to cleans incompatible characters from VIRSorter
    or VIRSorter2 affi headers.

    :param virsorter_in:

        Path to a affi tab file provide for DRAM by VIRSorter or VIRSorter2.
    :returns:
        String of file contents with headers cleaned.
    :raises ValueError:
        If the input format does not match VIRSorter or VIRSorter2.
    """
    if not is_affi_tab_not_fasta(virsorter_in):
        raise ValueError("The input file format does not match virsorter affi contigs, please check "
                         "that the flags match the file type: '-v' for virsorter, and '-i' for "
                         "fasta.")
    new_lines = list()
    for line in open(virsorter_in).readlines():
        line = line.strip()
        if line.startswith('>'):
            line = line.lstrip('>')
            split_line = line.split('|')
            split_line[0] = split_line[0].replace(';', '_').replace('=', '_')
            new_lines.append('>%s' % '|'.join(split_line))
        else:
            split_line = line.split('|')
            split_line[0] = split_line[0].replace(';', '_').replace('=', '_')
            new_lines.append('|'.join(split_line))
    return '%s\n' % '\n'.join(new_lines)


def remove_bad_chars(input_fasta: str = None, input_virsorter_affi_contigs: str = None, output: str = None):
    """
    Determine and call the function to cleans incompatible characters from headers of virsorter
    affi or fasta file.

    :param input_fasta:
        Path to a fasta file provide for DRAM by VIRSorter or VIRSorter2
    :param input_virsorter_affi_contigs:
        Path to a affi tab file provide for DRAM by VIRSorter
        or VIRSorter2
    :param output:
        Desired path for the new file clean of incompatible characters
    :raises ValueError:
        Check format, and invalid arguments.
    """
    if input_fasta is None and input_virsorter_affi_contigs is None:
        raise ValueError('Must provide either input fasta or virsorter affi contigs')
    elif input_fasta is not None:
        if ';' in output or '=' in output:
            raise ValueError('Cannot have = or ; in output fasta to work with DRAM-v, please change output fasta name')
        new_seqs = remove_bad_chars_fasta(input_fasta)
        write_sequence((seq for seq in new_seqs), format='fasta', into=output)
    else:
        virsorter_out = remove_bad_chars_virsorter_affi_contigs(input_virsorter_affi_contigs)
        open(output, 'w').write(virsorter_out)


def get_virsorter_hits(virsorter_affi_contigs):
    raw_file = open(virsorter_affi_contigs).read()

    virsorter_rows = list()
    for entry in raw_file[1:].split('\n>'):
        split_entry = entry.strip().split('\n')
        entry_data = split_entry[0].split('|')
        entry_name = entry_data[0]
        if '=' in entry_name or ';' in entry_name:
            raise ValueError('FASTA headers must not have = or ; before the first space (%s). To run DRAM-v you must '
                             'rerun VIRSorter with = and ; removed from the headers or run DRAM-v.py '
                             'remove_bad_characters and then rerun DRAM-v' % entry_name)
        entry_genes = split_entry[1:]
        entry_rows = [i.split('|') + [entry_name] for i in entry_genes]
        virsorter_rows += entry_rows
    virsorter_hits = pd.DataFrame(virsorter_rows, columns=VIRSORTER_COLUMN_NAMES).set_index('gene_name')
    virsorter_hits.index = [i.replace('(', '_').replace(')', '_') for i in virsorter_hits.index]
    virsorter_hits['name'] = [i.replace('(', '_').replace(')', '_') for i in virsorter_hits['name']]
    return virsorter_hits


def get_overlap(row1, row2):
    overlap_start = np.max((row1['start_position'], row2['start_position']))
    overlap_end = np.min((row1['end_position'], row2['end_position']))
    overlap_len = overlap_end - overlap_start
    if overlap_len <= 0:
        return 0, 0
    else:
        return overlap_len / (row1['end_position'] - row1['start_position']), \
               overlap_len / (row2['end_position'] - row2['start_position'])


def is_transposon(pfam_hits):
    if pd.isna(pfam_hits):
        return False
    else:
        pfams = {i[1:-1].split('.')[0] for i in re.findall(r'\[PF\d\d\d\d\d.\d*\]', pfam_hits)}
        return len(pfams & TRANSPOSON_PFAMS) > 0


def get_gene_order(dram_genes, virsorter_genes, min_overlap=.70):
    dram_genes['start_position'] = dram_genes['start_position'].astype(int)
    dram_genes['end_position'] = dram_genes['end_position'].astype(int)
    dram_genes = dram_genes.sort_values('start_position')
    dram_gene_number = 0

    virsorter_genes['start_position'] = virsorter_genes['start_position'].astype(int)
    virsorter_genes['end_position'] = virsorter_genes['end_position'].astype(int)
    virsorter_genes = virsorter_genes.sort_values('start_position')
    virsorter_gene_number = 0

    merged_genes_rows = list()
    while (dram_gene_number < dram_genes.shape[0]) and (virsorter_gene_number < virsorter_genes.shape[0]):
        dram_row = dram_genes.iloc[dram_gene_number]
        dram_gene = dram_genes.index[dram_gene_number]
        virsorter_row = virsorter_genes.iloc[virsorter_gene_number]
        virsorter_gene = virsorter_genes.index[virsorter_gene_number]
        virsorter_gene_category = virsorter_row['viral_protein_cluster_category']
        # end if we are at the end of either list
        forward_overlap, reverse_overlap = get_overlap(dram_row, virsorter_row)
        if (forward_overlap > min_overlap) and (reverse_overlap > min_overlap):
            merged_genes_rows.append((dram_gene, virsorter_gene, virsorter_gene_category))
            dram_gene_number += 1
            virsorter_gene_number += 1
        else:  # if not significantly overlap then determine first gene
            if dram_row['start_position'] < virsorter_row['start_position']:
                merged_genes_rows.append((dram_gene, None, None))
                dram_gene_number += 1
            elif dram_row['start_position'] > virsorter_row['start_position']:
                merged_genes_rows.append((None, virsorter_gene, virsorter_gene_category))
                virsorter_gene_number += 1
            else:  # if the same start positions then shorter gene is the first
                if dram_row['end_position'] < virsorter_row['end_position']:
                    merged_genes_rows.append((dram_gene, None, None))
                    dram_gene_number += 1
                elif dram_row['end_position'] > virsorter_row['end_position']:
                    merged_genes_rows.append((None, virsorter_gene,
                                              virsorter_gene_category))
                    virsorter_gene_number += 1
                else:  # genes have the same start and end position so this shouldn't happen
                    raise ValueError('This should never happen: %s, %s, %s, %s' % (dram_row.start_position,
                                                                                   dram_row.end_position,
                                                                                   virsorter_row.start_position,
                                                                                   virsorter_row.end_position))
    # clean up and add extras
    # if at end of both then just end
    if (dram_gene_number == dram_genes.shape[0]) and (virsorter_gene_number == virsorter_genes.shape[0]):
        pass
    # if not at end of dram annotations then add the rest
    elif dram_gene_number != dram_genes.shape[0]:
        for i in range(dram_gene_number, dram_genes.shape[0]):
            dram_gene = dram_genes.index[i]
            merged_genes_rows.append((dram_gene, None, None))
    # if not at end of virsorter genes then add the rest
    elif virsorter_gene_number != virsorter_genes.shape[0]:
        for i in range(virsorter_gene_number, virsorter_genes.shape[0]):
            virsorter_row = virsorter_genes.iloc[i]
            merged_genes_rows.append((None, virsorter_genes.index[i],
                                      virsorter_row['viral_protein_cluster_category']))
    # how can we be not at the end of either?
    else:
        return ValueError('How is this even possible?')
    return merged_genes_rows


def calculate_auxiliary_scores(gene_order):
    # right now saying that flanked means between here and the end of the contig
    gene_auxiliary_score_dict = dict()
    for i, (dram_gene, virsorter_gene, virsorter_category) in enumerate(gene_order):
        if dram_gene is not None:
            left_categories = [left_virsorter_category for left_dram_gene, left_viral_gene, left_virsorter_category
                               in gene_order[:i] if left_viral_gene is not None]
            hallmark_left = len(set(left_categories) & set(VIRSORTER_HALLMARK_GENE_CATEGORIES)) > 0
            viral_like_left = len(set(left_categories) & set(VIRSORTER_VIRAL_LIKE_GENE_CATEGORIES)) > 0

            right_categories = [right_virsorter_category for right_dram_gene, right_viral_gene, right_virsorter_category
                                in gene_order[i+1:] if right_viral_gene is not None]
            hallmark_right = len(set(right_categories) & set(VIRSORTER_HALLMARK_GENE_CATEGORIES)) > 0
            viral_like_right = len(set(right_categories) & set(VIRSORTER_VIRAL_LIKE_GENE_CATEGORIES)) > 0
            # if end of contig then 5:
            if i == 0 or i == (len(gene_order)-1):
                auxiliary_score = 5
            elif hallmark_left and hallmark_right:  # hallmark on both sides then cat 1
                auxiliary_score = 1
            # hallmark on one side and viral like on other then cat 2
            elif (hallmark_left and viral_like_right) or (viral_like_left and hallmark_right):
                auxiliary_score = 2
            # viral like on both side then cat 3
            elif viral_like_left and viral_like_right:
                auxiliary_score = 3
            # not end of contig and hallmark or viral like on other side
            elif (hallmark_left or viral_like_left) or (hallmark_right or viral_like_right):
                auxiliary_score = 4
            # if gene is the only virsorter viral like or hallmark on contig then make it a 4
            elif (virsorter_category in VIRSORTER_HALLMARK_GENE_CATEGORIES) or \
                 (virsorter_category in VIRSORTER_VIRAL_LIKE_GENE_CATEGORIES):
                auxiliary_score = 4
            else:  # if no viral like or hallmark genes then score is 5
                auxiliary_score = 5
            gene_auxiliary_score_dict[dram_gene] = auxiliary_score
    return gene_auxiliary_score_dict


def get_metabolic_flags(annotations, metabolic_genes, amgs, verified_amgs, scaffold_length_dict,
                        logger, length_from_end=5000):
    flag_dict = dict()
    metabolic_genes = set(metabolic_genes)
    for scaffold, scaffold_annotations in annotations.groupby('scaffold'):
        for gene, row in scaffold_annotations.iterrows():
            # set up
            flags = ''
            gene_annotations = set(get_ids_from_annotations_all(pd.DataFrame(row).transpose()).keys())  # TODO: Fix
            # is viral
            if 'vogdb_categories' in row.index and not pd.isna(row['vogdb_categories']):
                if len({'Xr', 'Xs'} & set(row['vogdb_categories'].split(';'))) > 0:
                    flags += 'V'
            # is metabolic
            if len(metabolic_genes & gene_annotations) > 0:
                flags += 'M'
            # is this a reported AMG reported
            if len(gene_annotations & set(amgs)) > 0:
                if 'M' not in flags:
                    flags += 'M'
                flags += 'K'
            # is this a experimentally verified amg
            if len(gene_annotations & set(verified_amgs)) > 0:
                flags += 'E'
            # is this gene a normal viral cell host entry gene
            if len(gene_annotations & CELL_ENTRY_CAZYS) > 0:
                flags += 'A'
            # is gene a normal virus peptidase
            if len(gene_annotations & VIRAL_PEPTIDASES_MEROPS) > 0:
                flags += 'P'
            # if there is a transposon in the contig
            if scaffold_annotations['is_transposon'].any():
                flags += 'T'
            # within 5 kb of end of contig
            if (int(row['start_position']) < length_from_end) or \
               (int(row['end_position']) > (scaffold_length_dict[row['scaffold']] - length_from_end)):
                flags += 'F'
            # if is_j:
            #     flags += 'J'
            flag_dict[gene] = flags
        # get 3 metabolic genes in a row flag
        for i in range(len(scaffold_annotations)):  # this needs to be fixed. Will only give B to middle of 3 genes.
            if 0 < i < (len(scaffold_annotations) - 1):
                gene = scaffold_annotations.index[i]
                gene_flags = flag_dict[gene]
                previous_gene = scaffold_annotations.index[i - 1]
                previous_gene_flags = flag_dict[previous_gene]
                next_gene = scaffold_annotations.index[i + 1]
                next_gene_flags = flag_dict[next_gene]
                if 'M' in previous_gene_flags and 'M' in gene_flags and 'M' in next_gene_flags:
                    if 'B' not in flag_dict[previous_gene]:
                        flag_dict[previous_gene] += 'B'
                    if 'B' not in flag_dict[gene]:
                        flag_dict[gene] += 'B'
                    if 'B' not in flag_dict[next_gene]:
                        flag_dict[next_gene] += 'B'
    return flag_dict


def get_amg_ids(amg_frame):
    ko_amgs = {j.strip() for i in amg_frame['KO'].dropna() for j in i.strip().split(';')}
    ec_amgs = {j.strip() for i in amg_frame['EC'].dropna() for j in i.strip().split(';')}
    pfam_amgs = {j.strip() for i in amg_frame['PFAM'].dropna() for j in i.strip().split(';')}
    # Put this back in once python 3.10 is more popular
    #return ko_amgs | ec_amgs | pfam_amgs
    return ko_amgs.union(ec_amgs, pfam_amgs)


def get_virsorter2_affi_contigs_name(scaffold):
    virsorter_scaffold_name = scaffold.replace('|', '_')
    return virsorter_scaffold_name


def get_virsorter_original_affi_contigs_name(scaffold):
    prophage_match = re.search(r'_gene_\d*_gene_\d*-\d*-\d*-cat_[123456]$', scaffold)
    plain_match = re.search(r'-cat_[123456]$', scaffold)
    if prophage_match is not None:
        virsorter_scaffold_name = scaffold[:prophage_match.start()]
    elif plain_match is not None:
        virsorter_scaffold_name = scaffold[:plain_match.start()]
    else:
        raise ValueError("Format is not VIRSorter2, and can't find VIRSorter "
                         "endings on fasta header: %s" % scaffold)
    return virsorter_scaffold_name


def get_virsorter_affi_contigs_name(scaffold):
    """
    :param scaffold: A raw string from the fasta header
    :returns: A string returned from the appropriate parser
    """
    virsorter2_match = re.search(r'\|\|', scaffold)
    if virsorter2_match is not None:
        return get_virsorter2_affi_contigs_name(scaffold)
    return get_virsorter_original_affi_contigs_name(scaffold)


def add_dramv_scores_and_flags(annotations, db_handler, logger, virsorter_hits=None, input_fasta=None):
    # setting up scoring viral genes
    amg_database_frame = pd.read_csv(db_handler.config['dram_sheets']['amg_database'], sep='\t')
    genome_summary_form = pd.read_csv(db_handler.config['dram_sheets']['genome_summary_form'], sep='\t', index_col=0)
    genome_summary_form = genome_summary_form.loc[genome_summary_form.potential_amg]

    # add auxiliary score
    if virsorter_hits is not None:
        gene_virsorter_category_dict = dict()
        gene_auxiliary_score_dict = dict()
        for scaffold, dram_frame in annotations.groupby('scaffold'):
            virsorter_scaffold_name = get_virsorter_affi_contigs_name(scaffold)
            virsorter_frame = virsorter_hits.loc[virsorter_hits.name == virsorter_scaffold_name]
            gene_order = get_gene_order(dram_frame, virsorter_frame)
            gene_virsorter_category_dict.update({dram_gene: virsorter_category for dram_gene, _, virsorter_category in
                                                 gene_order if dram_gene is not None})
            gene_auxiliary_score_dict.update(calculate_auxiliary_scores(gene_order))
        annotations['virsorter_category'] = pd.Series(gene_virsorter_category_dict)
        annotations['auxiliary_score'] = pd.Series(gene_auxiliary_score_dict)

    # get metabolic flags
    scaffold_length_dict = {seq.metadata['id']: len(seq) for seq in read_sequence(input_fasta, format='fasta')}
    metabolic_genes = set(genome_summary_form.index)
    if 'pfam_hits' in annotations:
        annotations['is_transposon'] = [is_transposon(i) for i in annotations['pfam_hits']]
    else:
        annotations['is_transposon'] = False

    amgs = get_amg_ids(amg_database_frame)
    verified_amgs = get_amg_ids(amg_database_frame.loc[amg_database_frame.verified])
    annotations['amg_flags'] = pd.Series(get_metabolic_flags(annotations, metabolic_genes, amgs,
                                                             verified_amgs, scaffold_length_dict,
                                                             logger))

    # downgrade B flag auxiliary scores
    if virsorter_hits is not None:
        annotations['auxiliary_score'] = pd.Series({gene: (4 if 'B' in row['amg_flags'] and row['auxiliary_score'] < 4
                                                           else row['auxiliary_score'])
                                                    for gene, row in annotations.iterrows()})

    return annotations

def annotate_vgfs(input_fasta, virsorter_affi_contigs=None, output_dir='.', min_contig_size=2500, split_contigs=False,
                  prodigal_mode='meta', trans_table='11', bit_score_threshold=60, rbh_bit_score_threshold=350,
                  custom_db_name=(), custom_fasta_loc=(), custom_hmm_loc=(), custom_hmm_cutoffs_loc=(),
                  custom_hmm_name=(), use_uniref=False, kofam_use_dbcan2_thresholds=False, skip_trnascan=False,
                  keep_tmp_dir=True, low_mem_mode=False, threads=10, verbose=True, config_loc:str=None ):
    mkdir(output_dir)
    log_file_path = path.join(output_dir, "annotate.log")
    logger = logging.getLogger('annotation_log')
    setup_logger(logger, log_file_path)

    # check inputs
    prodigal_modes = ['train', 'meta', 'single']
    if prodigal_mode not in prodigal_modes:
        raise ValueError('Prodigal mode must be one of %s.' % ', '.join(prodigal_modes))
    elif prodigal_mode in ['normal', 'single']:
        warnings.warn('When running prodigal in single mode your bins must have long contigs (average length >3 Kbp), '
                      'be long enough (total length > 500 Kbp) and have very low contamination in order for prodigal '
                      'training to work well.')

    # get database locations
    db_handler = DatabaseHandler(logger, config_loc=config_loc)
    db_handler.filter_db_locs(low_mem_mode=low_mem_mode,
                              use_uniref=use_uniref,
                              use_vogdb=True,
                              master_list=VMAG_DBS_TO_ANNOTATE)
    if virsorter_affi_contigs is not None:
        virsorter_hits = get_virsorter_hits(virsorter_affi_contigs)
    else:
        virsorter_hits = None

    # Check fastas
    perform_fasta_checks([input_fasta], logger)

    logger.info(f"Starting the Viral Annotation of Genes with database configuration: \n {db_handler.get_settings_str()}")

    if split_contigs:
        # split sequences into separate fastas
        contig_dir = path.join(output_dir, 'vMAGs')
        mkdir(contig_dir)

        contig_locs = list()
        for seq in read_sequence(input_fasta, format='fasta'):
            if len(seq) >= min_contig_size:
                if '=' in seq.metadata['id'] or ';' in seq.metadata['id']:
                    raise ValueError('FASTA headers must not have = or ; before the first space (%s). To run DRAM-v '
                                     'you must rerun VIRSorter with = and ; removed from the headers or run DRAM-v.py '
                                     'remove_bad_characters and then rerun DRAM-v' % seq.metadata['id'])
                if virsorter_hits is not None:
                    if get_virsorter_affi_contigs_name(seq.metadata['id']) not in virsorter_hits['name'].values:
                        raise ValueError("No virsorter calls found in %s for scaffold %s from input fasta" %
                                         (virsorter_affi_contigs, seq.metadata['id']))
                contig_loc = path.join(contig_dir, '%s.fasta' % seq.metadata['id'])
                write_sequence((i for i in [seq]), format='fasta', into=contig_loc)
                contig_locs.append(contig_loc)
    else:
        contig_locs = [input_fasta]

    # annotate vMAGs
    rename_bins = False
    annotations = annotate_fastas(
        fasta_locs=contig_locs,
        output_dir=output_dir, 
        db_handler=db_handler,
        logger=logger,
        min_contig_size=min_contig_size,
        prodigal_mode=prodigal_mode,
        trans_table=trans_table,
        bit_score_threshold=bit_score_threshold,
        rbh_bit_score_threshold=rbh_bit_score_threshold,
        custom_db_name=custom_db_name,
        custom_fasta_loc=custom_fasta_loc,
        custom_hmm_name=custom_hmm_name,
        custom_hmm_loc=custom_hmm_loc,
        custom_hmm_cutoffs_loc=custom_hmm_cutoffs_loc,
        kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
        skip_trnascan=skip_trnascan,
        rename_bins=rename_bins,
        keep_tmp_dir=keep_tmp_dir,
        threads=threads,
        verbose=verbose)
    logging.info('Annotations complete, assigning auxiliary scores and flags')

    annotations = add_dramv_scores_and_flags(annotations, db_handler, logger, virsorter_hits, input_fasta)

    # write annotations
    annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')

    logging.info("Completed annotations")

