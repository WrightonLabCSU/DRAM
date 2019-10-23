from datetime import datetime
from os import path, mkdir
from shutil import rmtree, copy2
import re

import pandas as pd
import numpy as np
from skbio.io import read as read_sequence

from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.annotate_bins import process_custom_dbs, annotate_fasta
from mag_annotator.utils import get_database_locs, get_ids_from_annotation

VIRSORTER_COLUMN_NAMES = ['gene_name', 'start_position', 'end_position', 'length', 'strandedness',
                          'viral_protein_cluster_hit', 'viral_protein_cluster_hit_score',
                          'viral_protein_cluster_hit_evalue', 'viral_protein_cluster_category', 'pfam_hit',
                          'pfam_hit_score', 'pfam_hit_evalue', 'name']

VIRSORTER_HALLMARK_GENE_CATEGORIES = {'0', '3'}

VIRSORTER_VIRAL_LIKE_GENE_CATEGORIES = {'1', '2', '4'}

TRANSPOSON_PFAMS = {'PF01609', 'PF00872', 'PF01610', 'PF01527', 'PF02371', 'PF01710', 'PF01385', 'PF01548', 'PF01526',
                    'PF01797', 'PF02899', 'PF05717', 'PF07592', 'PF03050', 'PF04754', 'PF04986', 'PF03400'}

CELL_ENTRY_CAZYS = {'CBM50', 'GH102', 'GH103', 'GH104', 'GH108', 'GH18', 'GH19', 'GH22', 'GH23', 'GH24', 'GH25', 'GH73',
                    'PL9', 'CBM12', 'CBM14', 'CBM18', 'CBM19'}


def get_virsorter_hits(virsorter_affi_contigs):
    raw_file = open(virsorter_affi_contigs).read()

    virsorter_rows = list()
    for entry in raw_file.split('>')[1:]:
        split_entry = entry.strip().split('\n')
        entry_data = split_entry[0].split('|')
        entry_name = entry_data[0]
        entry_genes = split_entry[1:]
        entry_rows = [i.split('|') + [entry_name] for i in entry_genes]
        virsorter_rows += entry_rows
    return pd.DataFrame(virsorter_rows, columns=VIRSORTER_COLUMN_NAMES).set_index('gene_name')


def get_overlap(row1, row2):
    overlap_start = np.max((row1['start_position'], row2['start_position']))
    overlap_end = np.min((row1['end_position'], row2['end_position']))
    overlap_len = overlap_end - overlap_start
    if overlap_len <= 0:
        return 0, 0
    else:
        return overlap_len / (row1.end_position - row1.start_position), \
               overlap_len / (row2.end_position - row2.start_position)


def get_next_number_name_row(i, frame):
    i += 1
    if i == frame.shape[0]:
        raise StopIteration()
    row = frame.iloc[i]
    row_name = row.name
    return i, row_name, row


def is_transposon(pfam_hits):
    if pd.isna(pfam_hits):
        return False
    else:
        pfams = {i[1:-1].split('.')[0] for i in re.findall('\[PF\d\d\d\d\d.\d*\]', pfam_hits)}
        return len(pfams & TRANSPOSON_PFAMS) > 0


def get_gene_order(dram_genes, virsorter_genes, min_overlap=.70):
    dram_gene_frame_list = list()
    dram_gene_frame_list.append(dram_genes['start_position'].astype(int))
    dram_gene_frame_list.append(dram_genes['end_position'].astype(int))
    dram_gene_frame = pd.DataFrame(dram_gene_frame_list).transpose()
    dram_gene_frame = dram_gene_frame.sort_values('start_position')
    dram_gene_number, dram_gene, dram_row = get_next_number_name_row(-1, dram_gene_frame)

    virsorter_gene_frame_list = list()
    virsorter_gene_frame_list.append(virsorter_genes['start_position'].astype(int))
    virsorter_gene_frame_list.append(virsorter_genes['end_position'].astype(int))
    virsorter_gene_frame_list.append(virsorter_genes['viral_protein_cluster_category'])
    virsorter_gene_frame = pd.DataFrame(virsorter_gene_frame_list).transpose()
    virsorter_gene_frame = virsorter_gene_frame.sort_values('start_position')
    virsorter_gene_number, virsorter_gene, virsorter_row = get_next_number_name_row(-1, virsorter_gene_frame)

    merged_genes_rows = list()
    while True:
        try:
            # end if we are at the end of either list
            forward_overlap, reverse_overlap = get_overlap(dram_row, virsorter_row)
            if (forward_overlap > min_overlap) and (reverse_overlap > min_overlap):
                merged_genes_rows.append((dram_gene, virsorter_gene, virsorter_row['viral_protein_cluster_category']))
                dram_gene_number, dram_gene, dram_row = get_next_number_name_row(dram_gene_number, dram_gene_frame)
                virsorter_gene_number, virsorter_gene, virsorter_row = get_next_number_name_row(virsorter_gene_number,
                                                                                                virsorter_gene_frame)
            else:  # if not significantly overlap then determine first gene
                if dram_row['start_position'] < virsorter_row['start_position']:
                    merged_genes_rows.append((dram_gene, None, None))
                    dram_gene_number, dram_gene, dram_row = get_next_number_name_row(dram_gene_number, dram_gene_frame)
                elif dram_row['start_position'] > virsorter_row['start_position']:
                    merged_genes_rows.append((None, virsorter_gene, virsorter_row['viral_protein_cluster_category']))
                    virsorter_gene_number, virsorter_gene, virsorter_row = \
                        get_next_number_name_row(virsorter_gene_number, virsorter_gene_frame)
                else:  # if the same start positions then shorter gene is the first
                    if dram_row['end_position'] < virsorter_row['end_position']:
                        merged_genes_rows.append((dram_gene, None, None))
                        dram_gene_number, dram_gene, dram_row = get_next_number_name_row(dram_gene_number,
                                                                                         dram_gene_frame)
                    elif dram_row['end_position'] > virsorter_row['end_position']:
                        merged_genes_rows.append((None, virsorter_gene,
                                                  virsorter_row['viral_protein_cluster_category']))
                        virsorter_gene_number, virsorter_gene, virsorter_row = \
                            get_next_number_name_row(virsorter_gene_number, virsorter_gene_frame)
                    else:  # genes have the same start and end position so this shouldn't happen
                        raise ValueError('This should never happen: %s, %s, %s, %s' % (dram_row.start_position,
                                                                                       dram_row.end_position,
                                                                                       virsorter_row.start_position,
                                                                                       virsorter_row.end_position))
        except StopIteration:
            break
    # clean up and add extras
    # if at end of both then just end
    if (dram_gene_number == dram_gene_frame.shape[0]-1) and (virsorter_gene_number == virsorter_gene_frame.shape[0]-1):
        pass
    # if not at end of dram annotations then add the rest
    elif dram_gene_number != dram_genes.shape[0]-1:
        for i in range(dram_gene_number, dram_gene_frame.shape[0]):
            merged_genes_rows.append((dram_gene_frame.index[i], None, None))
    # if not at end of virsorter genes then add the rest
    elif virsorter_gene_number != virsorter_gene_frame.shape[0]-1:
        for i in range(virsorter_gene_number, virsorter_gene_frame.shape[0]):
            merged_genes_rows.append((None, virsorter_gene_frame.index[i],
                                      virsorter_gene_frame.viral_protein_cluster_category[i]))
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
                                in gene_order[i + 1:] if right_viral_gene is not None]
            hallmark_right = len(set(right_categories) & set(VIRSORTER_HALLMARK_GENE_CATEGORIES)) > 0
            viral_like_right = len(set(right_categories) & set(VIRSORTER_VIRAL_LIKE_GENE_CATEGORIES)) > 0
            if hallmark_left and hallmark_right:  # hallmark on both sides then cat 1
                auxiliary_score = 1
            # hallmark on one side and viral like on other then cat 2
            elif (hallmark_left or viral_like_left) and (hallmark_right or viral_like_right):
                auxiliary_score = 2
            # viral like on both side then cat 3
            elif viral_like_left and viral_like_right:
                auxiliary_score = 3
            # not end of contig and hallmark or viral like on other side
            elif ((hallmark_left or viral_like_left) and len(right_categories) > 0) or \
                    (len(left_categories) > 0 and (hallmark_right or viral_like_right)):
                auxiliary_score = 4
            else:  # if gene is at end of contig or no viral like or hallmark genes then score is 5
                auxiliary_score = 5
            print(i, auxiliary_score, hallmark_left, viral_like_left, hallmark_right, viral_like_right)
            gene_auxiliary_score_dict[dram_gene] = auxiliary_score
    return gene_auxiliary_score_dict


def get_metabolic_flags(annotations, metabolic_genes, amgs, verified_amgs, scaffold_length_dict,
                        length_from_end=5000):
    flag_dict = dict()
    metabolic_genes = set(metabolic_genes)
    for scaffold, scaffold_annotations in annotations.groupby('scaffold'):
        for gene, row in scaffold_annotations.iterrows():
            # set up
            flags = ''
            gene_annotations = set(get_ids_from_annotation(pd.DataFrame(row).transpose()).keys())
            # is viral
            if not pd.isna(row['vogdb_categories']):
                if len({'Xr', 'Xs'} & set(row['vogdb_categories'].split(';'))) > 0:
                    flags += 'V'
            # is metabolic
            if len(metabolic_genes & gene_annotations) > 0:
                flags += 'M'
            # is this a reported AMG reported
            if len(gene_annotations & set(amgs)) > 0:
                flags += 'K'
            # is this a experimentally verified amg
            if len(gene_annotations & set(verified_amgs)) > 0:
                flags += 'E'
            # is this gene a normal viral cell host entry gene
            if len(gene_annotations & CELL_ENTRY_CAZYS) > 0:
                flags += 'A'
            # if there is a transposon in the contig
            if scaffold_annotations['is_transposon'].any():
                flags += 'T'
            # within 5 kb of end of contig
            if (int(row.start_position) < length_from_end) or \
                    (int(row.end_position) > (scaffold_length_dict[row.scaffold] - length_from_end)):
                flags += 'F'
            flag_dict[gene] = flags
        # get 3 metabolic genes in a row flag
        for i in range(len(scaffold_annotations)):
            if 0 < i < (len(scaffold_annotations) - 1):
                gene = scaffold_annotations.index[i]
                gene_flags = flag_dict[gene]
                previous_gene = scaffold_annotations.index[i - 1]
                previous_gene_flags = flag_dict[previous_gene]
                next_gene = scaffold_annotations.index[i + 1]
                next_gene_flags = flag_dict[next_gene]
                if 'M' in previous_gene_flags and 'M' in gene_flags and 'M' in next_gene_flags:
                    flag_dict[gene] += 'B'
    return flag_dict


def get_amg_ids(amg_frame):
    ko_amgs = set([j.strip() for i in amg_frame['KO'].dropna() for j in i.strip().split(';')])
    ec_amgs = set([j.strip() for i in amg_frame['EC'].dropna() for j in i.strip().split(';')])
    pfam_amgs = set([j.strip() for i in amg_frame['PFAM'].dropna() for j in i.strip().split(';')])
    return ko_amgs | ec_amgs | pfam_amgs


def get_virsorter_affi_contigs_name(scaffold):
    prophage_match = re.search(r'_gene_\d*_gene_\d*-\d*-\d*-cat_[1,2,3,4,5,6]$', scaffold)
    plain_match = re.search(r'-cat_[1,2,3,4,5,6]$', scaffold)
    if prophage_match is not None:
        virsorter_scaffold_name = scaffold[:prophage_match.start()]
    elif plain_match is not None:
        virsorter_scaffold_name = scaffold[:plain_match.start()]
    else:
        raise ValueError("Can't find VIRSorter endings on fasta header: %s" % scaffold)
    return virsorter_scaffold_name


def annotate_vgfs(input_fasta, virsorter_affi_contigs, output_dir='.', min_contig_size=5000, bit_score_threshold=60,
                  rbh_bit_score_threshold=350, custom_db_name=(), custom_fasta_loc=(), skip_uniref=True,
                  skip_trnascan=False, keep_tmp_dir=True, threads=10, verbose=True):
    # set up
    start_time = datetime.now()
    print('%s: Annotation started' % str(datetime.now()))

    mkdir(output_dir)
    tmp_dir = path.join(output_dir, 'working_dir')
    mkdir(tmp_dir)

    # get database locations
    db_locs = get_database_locs()
    db_handler = DatabaseHandler(db_locs['description_db'])
    custom_db_locs = process_custom_dbs(custom_fasta_loc, custom_db_name, tmp_dir, threads, verbose)
    print('%s: Retrieved database locations and descriptions' % (str(datetime.now() - start_time)))

    # iterate over list of fastas and annotate each individually
    # get name of file e.g. /home/shaffemi/my_genome.fa -> my_genome
    fasta_name = path.splitext(path.basename(input_fasta.strip('.gz')))[0]
    annotations = annotate_fasta(input_fasta, fasta_name, tmp_dir, db_locs, db_handler, min_contig_size, custom_db_locs,
                                 None, bit_score_threshold, rbh_bit_score_threshold, skip_uniref, skip_trnascan,
                                 start_time, threads, verbose)
    print('%s: Annotations complete, processing annotations' % str(datetime.now() - start_time))

    # setting up scoring viral genes
    amg_database_frame = pd.read_csv(db_locs['amg_database'], sep='\t')
    genome_summary_frame = pd.read_csv(db_locs['genome_summary_form'], sep='\t', index_col=0)
    virsorter_hits = get_virsorter_hits(virsorter_affi_contigs)

    # add auxiliary score
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
    metabolic_genes = set(genome_summary_frame.index)
    annotations['is_transposon'] = [is_transposon(i) for i in annotations['pfam_hits']]

    amgs = get_amg_ids(amg_database_frame)
    verified_amgs = get_amg_ids(amg_database_frame.loc[amg_database_frame.verified])
    annotations['amg_flags'] = pd.Series(get_metabolic_flags(annotations, metabolic_genes, amgs,
                                                             verified_amgs, scaffold_length_dict))

    # write annotations
    annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')

    # copy results files to output
    copy2(path.join(tmp_dir, 'genes.annotated.fna'), path.join(output_dir, 'genes.fna'))
    copy2(path.join(tmp_dir, 'genes.annotated.faa'), path.join(output_dir, 'genes.faa'))
    copy2(path.join(tmp_dir, 'scaffolds.annotated.fa'), path.join(output_dir, 'scaffolds.fna'))
    copy2(path.join(tmp_dir, 'genes.annotated.gff'), path.join(output_dir, 'genes.gff'))
    copy2(path.join(tmp_dir, 'trnas.tsv'), path.join(output_dir, 'trnas.tsv'))
    copy2(path.join(tmp_dir, 'rrnas.tsv'), path.join(output_dir, 'rrnas.tsv'))
    copy2(path.join(tmp_dir, '%s.gbk' % fasta_name), path.join(output_dir, 'scaffolds.gbk'))

    # clean up
    if not keep_tmp_dir:
        rmtree(tmp_dir)

    print("%s: Completed annotations" % str(datetime.now() - start_time))
