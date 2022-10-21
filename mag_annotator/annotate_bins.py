"""
Main control point for the annotation process
"""
import re
import io
import time
import logging
from glob import glob
from functools import partial
from datetime import datetime
from skbio.io import read as read_sequence
from skbio.io import write as write_sequence
from skbio import Sequence
from skbio.metadata import IntervalMetadata
from os import path, mkdir, stat
from shutil import rmtree, copy2
import pandas as pd

from mag_annotator.utils import run_process, make_mmseqs_db, merge_files, \
    multigrep, remove_suffix, setup_logger, run_hmmscan, sig_scores, get_sig_row, \
    generic_hmmscan_formater, get_reciprocal_best_hits, get_best_hits, BOUTFMT6_COLUMNS
from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.fasta_dup_name_test import fastas_dup_check

MAG_DBS_TO_ANNOTATE = ('kegg', 'kofam_hmm', 'kofam_ko_list', 'uniref', 'peptidase', 'pfam', 'dbcan', 'vogdb')

# TODO: Bind verbose to logging
# TODO Exceptions are not fully handled by logging


def filter_fasta(fasta_loc, min_len=5000, output_loc=None):
    """Removes sequences shorter than a set minimum from fasta files, outputs an object or to a file"""
    kept_seqs = (seq for seq in read_sequence(fasta_loc, format='fasta') if len(seq) >= min_len)
    if output_loc is None:
        return list(kept_seqs)
    else:
        write_sequence(kept_seqs, format='fasta', into=output_loc)


def run_prodigal(fasta_loc, output_dir, logger, mode='meta', trans_table='11', verbose=False):
    """Runs the prodigal gene caller on a given fasta file, outputs resulting files to given directory"""
    output_gff = path.join(output_dir, 'genes.gff')
    output_fna = path.join(output_dir, 'genes.fna')
    output_faa = path.join(output_dir, 'genes.faa')

    run_process(['prodigal', '-i', fasta_loc, '-p', mode, '-g', trans_table, '-f', 'gff', '-o', output_gff, '-a',
                 output_faa, '-d', output_fna], logger, verbose=verbose)
    return output_gff, output_fna, output_faa


def process_reciprocal_best_hits(forward_output_loc, reverse_output_loc, target_prefix='target'):
    """Process the forward and reverse best hits results to find reverse best hits
    Returns the query gene, target gene, if it was a reverse best hit, % identity, bit score and e-value
    """
    forward_hits = pd.read_csv(forward_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    forward_hits = forward_hits.set_index('qId')
    reverse_hits = pd.read_csv(reverse_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    reverse_hits = reverse_hits.set_index('qId')

    def check_hit(row:pd.Series):
        rbh = False
        if row.tId in reverse_hits.index:
            rbh = row.name == reverse_hits.loc[row.tId].tId
        return {'%s_hit' % target_prefix:      row.tId,
                '%s_RBH' % target_prefix:      rbh,
                '%s_identity' % target_prefix: row.seqIdentity,
                '%s_bitScore' % target_prefix: row.bitScore,
                '%s_eVal' % target_prefix:     row.eVal,
                'index':                       row.name
                }
    hits = forward_hits.apply(check_hit, axis=1, result_type='expand')
    # NOTE these lines may not be necessary
    hits.set_index('index', drop=True, inplace=True)
    hits.index.name = None
    return hits


def get_kegg_description(kegg_hits, header_dict):
    """Gets the KEGG IDs, and full KEGG hits from list of KEGG IDs for output in annotations"""
    gene_description = list()
    ko_list = list()
    for kegg_hit in kegg_hits.kegg_hit:
        header = header_dict[kegg_hit]
        gene_description.append(header)
        kos = re.findall(r'(K\d\d\d\d\d)', header)
        if len(kos) == 0:
            ko_list.append('')
        else:
            ko_list.append(','.join(kos))
    # TODO: change kegg_id to kegg_genes_id so that people get an error and not the wrong identifier
    new_df = pd.DataFrame([kegg_hits['kegg_hit'].values, ko_list, gene_description],
                          index=['kegg_genes_id', 'ko_id', 'kegg_hit'], columns=kegg_hits.index)
    return pd.concat([new_df.transpose(), kegg_hits.drop('kegg_hit', axis=1)], axis=1, sort=False)


def get_uniref_description(uniref_hits, header_dict):
    """Gets UniRef ID's, taxonomy and full string from list of UniRef IDs for output in annotations"""
    gene_description = list()
    uniref_list = list()
    gene_taxonomy = list()
    for uniref_hit in uniref_hits.uniref_hit:
        header = header_dict[uniref_hit]
        gene_description.append(header)
        uniref_list.append(header[header.find('RepID=') + 6:])
        gene_taxonomy.append(re.search(r'Tax=(.*?) (\S*?)=', header).group(1))
    new_df = pd.DataFrame([uniref_list, gene_description, gene_taxonomy],
                          index=['uniref_id', 'uniref_hit', 'uniref_taxonomy'],
                          columns=uniref_hits.index)
    return pd.concat([new_df.transpose(), uniref_hits.drop('uniref_hit', axis=1)], axis=1, sort=False)


def get_basic_description(hits, header_dict, db_name='viral'):
    """Get viral gene full descriptions based on headers (text before first space)"""
    hit_list = list()
    description = list()
    [hit for hit in hits['%s_hit' % db_name]]
    for hit in hits['%s_hit' % db_name]:
        header = header_dict[hit]
        hit_list.append(hit)
        description.append(header)
    new_df = pd.DataFrame([hit_list, description],
                          index=['%s_id' % db_name, '%s_hit' % db_name],
                          columns=hits.index)
    return pd.concat([new_df.transpose(), hits.drop('%s_hit' % db_name, axis=1)], axis=1, sort=False)


def get_peptidase_description(peptidase_hits, header_dict):
    peptidase_list = list()
    peptidase_family = list()
    peptidase_descirption = list()
    for peptidase_hit in peptidase_hits.peptidase_hit:
        header = header_dict[peptidase_hit]
        peptidase_list.append(peptidase_hit)
        peptidase_family.append(re.search(r'#\w*.#', header).group()[1:-1])
        peptidase_descirption.append(header)
    new_df = pd.DataFrame([peptidase_list, peptidase_family, peptidase_descirption],
                          index=['peptidase_id', 'peptidase_family', 'peptidase_hit'], columns=peptidase_hits.index)
    return pd.concat([new_df.transpose(), peptidase_hits.drop('peptidase_hit', axis=1)], axis=1, sort=False)


def run_mmseqs_profile_search(query_db, pfam_profile, output_loc, logger, output_prefix='mmpro_results', db_handler=None,
                              threads=10, verbose=False):
    """Use mmseqs to run a search against pfam, currently keeping all hits and not doing any extra filtering"""
    tmp_dir = path.join(output_loc, 'tmp')
    output_db = path.join(output_loc, '%s.mmsdb' % output_prefix)
    run_process(['mmseqs', 'search', query_db, pfam_profile, output_db, tmp_dir, '-k', '5', '-s', '7', '--threads',
                 str(threads)], logger, verbose=verbose)
    output_loc = path.join(output_loc, '%s_output.b6' % output_prefix)
    run_process(['mmseqs', 'convertalis', query_db, pfam_profile, output_db, output_loc], logger, verbose=verbose)
    pfam_results = pd.read_csv(output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    if pfam_results.shape[0] > 0:
        pfam_dict = dict()
        if db_handler is not None:
            pfam_descriptions = db_handler.get_descriptions(set(pfam_results.tId), '%s_description' % output_prefix)
        else:
            pfam_descriptions = {}
        for gene, pfam_frame in pfam_results.groupby('qId'):
            if len(pfam_descriptions) < 1:
                pfam_dict[gene] = '; '.join(pfam_frame.tId)
            else:
                pfam_dict[gene] = '; '.join(['%s [%s]' % (pfam_descriptions[ascession], ascession)
                                             for ascession in pfam_frame.tId])
        return pd.DataFrame(pfam_dict, index=[f"{output_prefix}_hits"]).T
    else:
        return pd.DataFrame(columns=[f"{output_prefix}_hits"])


def find_best_dbcan_hit(genome:str, group:pd.DataFrame):
    group['perc_cov'] = group.apply(
        lambda x: (x['target_end'] - x['target_start']) / x['target_length'],
        axis=1)
    group.sort_values('perc_cov', inplace=True)
    group.columns
    group.sort_values('full_evalue', inplace=True)
    return group.iloc[0]['target_id']


def dbcan_hmmscan_formater(hits:pd.DataFrame,  db_name:str, db_handler=None):
    """
    format the ouput of the dbcan database.

    Note these changes
    introduce new column for best hit from cazy database- this will be the best hit above the already established threshold (0.35 coverage, e-18 evalue) - then for the distillate pull info from the best hit only
    introduce new column for corresponding EC number information from sub families (EC numbers are subfamily ECs)
    Make sure that ids and descriptions are separate (descriptions these are family based)
    :param hits:
    :param db_name:
    :param db_handler:
    :returns:
    """
    hits_sig = hits[hits.apply(partial(get_sig_row, evalue_lim=1e-18), axis=1)]
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    hit_groups = hits_sig.groupby('query_id')
    all_hits = hit_groups.apply(
        lambda x: '; '.join(x['target_id'].apply(lambda y:y[:-4]).unique())
    )
    hits_df = pd.DataFrame(all_hits)
    hits_df.columns = [f"{db_name}_ids"]
    def description_pull(x:str):
         id_list = [re.findall('^[A-Z]*[0-9]*', str(x))[0] for x in  x.split('; ')],
         id_list = [y for x in  id_list for y in x if len(x) > 0]
         description_list = db_handler.get_descriptions(id_list, 'dbcan_description').values()
         description_str = '; '.join(description_list)
         return description_str
    if db_handler is not None:
        hits_df[f"{db_name}_hits"] = hits_df[f"{db_name}_ids"].apply(description_pull)
        hits_df[f"{db_name}_subfam_ec"] = hits_df[f"{db_name}_ids"].apply(
            lambda x:'; '.join(
                db_handler.get_descriptions(
                    x.split('; '),
                    'dbcan_description',
                    description_name='ec'
                ).values()))
    hits_df[f'{db_name}_best_hit'] = [find_best_dbcan_hit(*i) for i in hit_groups]
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns
    return hits_df


def kofam_hmmscan_formater(hits:pd.DataFrame, hmm_info_path:str=None, use_dbcan2_thresholds:bool=False,
                           top_hit:bool=True):
    hmm_info = pd.read_csv(hmm_info_path, sep='\t', index_col=0)
    if use_dbcan2_thresholds:
        hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hits_sig = sig_scores(hits, hmm_info)
    # if there are any significant results then parse to dataframe
    if len(hits_sig) == 0:
        return pd.DataFrame()
    kegg_dict = dict()
    for gene, frame in hits_sig.groupby('query_id'):
        # TODO: take top hit for full length genes and all hits for domains?
        # TODO: if top hit then give all e-value and bitscore info
        if top_hit:
            best_hit = frame[frame.full_evalue == frame.full_evalue.min()]
            ko_id = best_hit['target_id'].iloc[0]
            kegg_dict[gene] = [ko_id, hmm_info.loc[ko_id, 'definition']]
        else:
            kegg_dict[gene] = [','.join([i for i in frame.target_id]),
                               '; '.join([hmm_info.loc[i, 'definition'] for i in frame.target_id])]
    return pd.DataFrame(kegg_dict, index=['ko_id', 'kegg_hit']).transpose()



def vogdb_hmmscan_formater(hits:pd.DataFrame,  db_name:str, logger:logging.Logger,
                           db_handler=None):
    categories_col = f"{db_name}_categories"
    id_col = f"{db_name}_id"
    description_col = f"{db_name}_description"
    hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    if len(hits_sig) == 0:
        logger.warn("No significant hits for vog_db")
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame(columns=[categories_col, id_col, description_col])
    # Get the best hits
    hits_best = hits_sig.sort_values('full_evalue').drop_duplicates(subset=["query_id"])
    if db_handler is None:
        hits_df = hits_best[['target_id', 'query_id']].copy()
    else:
        # get_descriptions
        desc_col = f"{db_name}_hits"
        descriptions = pd.DataFrame(
            db_handler.get_descriptions(hits_best['target_id'].unique(), description_col),
            index=[desc_col]).T
        categories = descriptions[desc_col].apply(lambda x: x.split('; ')[-1])
        descriptions[categories_col] = categories.apply(
            lambda x: ';'.join(set([x[i:i + 2] for i in range(0, len(x), 2)])))
        descriptions['target_id'] = descriptions.index
        hits_df = pd.merge(hits_best[['query_id', 'target_id']], descriptions, on=f'target_id')
    hits_df.set_index('query_id', inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    hits_df.rename(columns={'target_id': id_col}, inplace=True)
    return hits_df


def get_gene_data(fasta_loc):
    """Take the prodigal gene headers and get the scaffold that it came from
    Based on idba_ud 'scaffold_#' scaffold names with gene name after
    """
    df_dict = dict()
    for seq in read_sequence(fasta_loc, format='fasta'):
        split_label = seq.metadata['id'].split('_')
        scaffold = '_'.join(split_label[:-1])
        gene_position = split_label[-1]
        start_position, end_position, strandedness = seq.metadata['description'].split('#')[1:4]
        df_dict[seq.metadata['id']] = [scaffold, int(gene_position), int(start_position), int(end_position),
                                       int(strandedness)]
    return pd.DataFrame.from_dict(df_dict, orient='index', columns=['scaffold', 'gene_position', 'start_position',
                                                                    'end_position', 'strandedness'])


def get_unannotated(fasta_loc, annotations):
    """Get the genes from the fasta which did not get any annotations"""
    return [seq.metadata['id'] for seq in read_sequence(fasta_loc, format='fasta')
            if seq.metadata['id'] not in annotations]


def assign_grades(annotations):
    """Grade genes based on reverse best hits to KEGG, UniRef and Pfam"""
    grades = dict()
    for gene, row in annotations.iterrows():
        if row.get('kegg_RBH') is True:
            rank = 'A'
        elif row.get('uniref_RBH') is True:
            rank = 'B'
        elif not pd.isna(row.get('kegg_hit')) or not pd.isna(row.get('uniref_hit')):
            rank = 'C'
        elif not pd.isna(row.get('pfam_hits')) or not pd.isna(row.get('cazy_hits'))\
                or not pd.isna(row.get('peptidase_hit')):
            rank = 'D'
        else:
            rank = 'E'
        grades[gene] = rank
    return pd.DataFrame(grades, index=['rank']).T


def generate_annotated_fasta(input_fasta, annotations, verbosity='short', name=None):
    """Generates fasta entries with added annotation information to the header of a fasta
    either add best annotation (based on grade) (verbosity = short) or all annotations (verbosity = long)
    """
    for seq in read_sequence(input_fasta, format='fasta'):
        annotation = annotations.loc[seq.metadata['id']]
        if 'rank' in annotations.columns and verbosity == 'short':
            annotation_str = 'rank: %s' % annotation['rank']
            if annotation['rank'] == 'A':
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'kegg')
            elif annotation['rank'] == 'B':
                annotation_str += '; %s (db=%s)' % (annotation.uniref_hit, 'uniref')
            elif annotation['rank'] == 'C':
                if 'kegg_hit' in annotation:
                    annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'kegg')
                if 'uniref_hit' in annotation:
                    annotation_str += '; %s (db=%s)' % (annotation.uniref_hit, 'uniref')
            elif annotation['rank'] == 'D':
                annotation_str += '; %s (db=%s)' % (annotation.pfam_hits, 'pfam')
            else:
                pass
        else:
            annotation_list = []
            if 'rank' in annotations.columns:
                annotation_list += ['rank: %s' % annotation['rank']]
            if 'kegg_hit' in annotations.columns:
                if not pd.isna(annotation.kegg_hit):
                    annotation_list += ['%s (db=%s)' % (annotation.kegg_hit, 'kegg')]
            if 'uniref_hit' in annotations.columns:
                if not pd.isna(annotation.uniref_hit):
                    annotation_list += ['%s (db=%s)' % (annotation.uniref_hit, 'uniref')]
            if 'pfam_hits' in annotations.columns:
                if not pd.isna(annotation.pfam_hits):
                    annotation_list += ['%s (db=%s)' % (annotation.pfam_hits, 'pfam')]
            if 'bin_taxonomy' in annotations.columns:
                if not pd.isna(annotation.bin_taxonomy):
                    annotation_list += [annotation.bin_taxonomy]
            annotation_str = '; '.join(annotation_list)
        if name is not None:
            seq.metadata['id'] = '%s_%s' % (name, seq.metadata['id'])
        seq.metadata['description'] = annotation_str
        yield seq


def create_annotated_fasta(input_fasta, annotations, output_fasta, verbosity='short', name=None):
    """For use with genes files, added annotations"""

    write_sequence( generate_annotated_fasta(input_fasta, annotations, verbosity, name), format='fasta', into=output_fasta)


def generate_renamed_fasta(input_fasta, prefix):
    """For use with scaffolds files, merges together bins with fasta name added as a prefix to the file"""
    for seq in read_sequence(input_fasta, format='fasta'):
        seq.metadata['id'] = '%s_%s' % (prefix, seq.metadata['id'])
        yield seq


def rename_fasta(input_fasta, output_fasta, prefix):
    """See above"""
    write_sequence(generate_renamed_fasta(input_fasta, prefix), format='fasta', into=output_fasta)


# TODO: add annotations that don't end with '_id'
def annotate_gff(input_gff, output_gff, annotations, prefix=None):
    """Go through a gff and add a prefix to the scaffold and gene number for all ID's"""
    f = open(input_gff)
    o = open(output_gff, 'w')
    for line in f:
        line = line.strip()
        if not line.startswith('#'):
            # replace id with new name
            old_scaffold = line.split('\t')[0]
            match = re.search(r'ID=\d*_\d*;', line)
            gene_number = match.group().split('_')[-1][:-1]
            old_gene_name = '%s_%s' % (old_scaffold, gene_number)
            if prefix is not None:  # add prefix to line name and gene name if given
                line = '%s_%s' % (prefix, line)
                gene_name = '%s_%s' % (prefix, old_gene_name)
            else:
                gene_name = old_gene_name
            line = re.sub(r'ID=\d*_\d*;', 'ID=%s;' % gene_name, line)
            # get annotations to add from annotations file and add to end of line
            annotations_to_add = {strip_endings(i, ['_id']): annotations.loc[old_gene_name, i]
                                  for i in annotations.columns if i.endswith('_id')}
            database_information = ['Dbxref="%s:%s"' % (key.strip(), value.strip().replace(':', '_'))
                                    for key, values in annotations_to_add.items()
                                    if not pd.isna(values)
                                    for value in values.split('; ')
                                    if value != '']
            if len(database_information) > 0:
                line += '%s;' % ';'.join(database_information)
        o.write('%s\n' % line)


def make_gbk_from_gff_and_fasta(gff_loc='genes.gff', fasta_loc='scaffolds.fna', faa_loc='genes.faa', output_gbk=None):
    # filter scaffolds out of fasta which are not in genes.gff
    gff_frame = pd.read_csv(gff_loc, sep='\t', comment='#', header=None)
    gff_scaffolds = set(gff_frame[0])
    capture_fasta = io.StringIO()
    write_sequence((i for i in read_sequence(fasta_loc, format='fasta') if i.metadata['id'] in gff_scaffolds),
                   into=capture_fasta, format='fasta')
    # scikit-bio can make a genbank for a concatenated gff and fasta file so make this first
    gff = open(gff_loc)
    aa_records = {i.metadata['id']: i for i in read_sequence(faa_loc, format='fasta')}
    # the gff with fasta format is the gff then ##FASTA then the fasta
    concat_gff = '%s\n##FASTA\n%s' % (gff.read(), capture_fasta.getvalue())

    # now we can only ready one record from the gff/fasta hybrid at a time
    # read a record at a time till end of the fasta
    genbank_records = ''
    # +1 because there is no \n> for the first line in the file and +1 because we are indexing starting at 1
    for i in range(1, capture_fasta.getvalue().count('\n>') + 1 + 1):
        seq = read_sequence(io.StringIO(concat_gff), format='gff3', into=Sequence, seq_num=i)
        seq.metadata['LOCUS'] = {'locus_name': seq.metadata['id'],
                                 'size': len(seq),
                                 'unit': 'bp',
                                 'mol_type': 'DNA',
                                 'shape': 'linear',
                                 'division': 'ENV',
                                 'date': time.strftime('%d-%b-%Y').upper()}
        seq_bounds = (seq.interval_metadata.lower_bound, seq.interval_metadata.upper_bound)
        for interval in seq.interval_metadata.query([seq_bounds]):
            interval.metadata['gene'] = interval.metadata['ID']
            if interval.metadata['ID'] in aa_records:
                interval.metadata['translation'] = aa_records[interval.metadata['ID']]
        # need to capture genbank output so we can combine into a multi-genbank file
        capture_print = io.StringIO()
        seq.write(capture_print, format='genbank')
        genbank_records += capture_print.getvalue()

    # if no output then return as string, otherwise write write to output location
    if output_gbk is None:
        return genbank_records
    else:
        open(output_gbk, 'w').write(genbank_records)


def get_dups(columns):
    keep = list()
    seen = list()
    for column in columns:
        if column in seen:
            keep.append(False)
        else:
            seen.append(column)
            keep.append(True)
    return keep


def run_trna_scan(fasta, tmp_dir, fasta_name, logger, threads=10, verbose=True):
    """Run tRNAscan-SE on scaffolds and create a table of tRNAs as a separate output"""
    raw_trnas = path.join(tmp_dir, 'raw_trnas.txt')
    run_process(['tRNAscan-SE', '-G', '-o', raw_trnas, '--thread', str(threads), fasta], logger, verbose=verbose)
    if path.isfile(raw_trnas) and stat(raw_trnas).st_size > 0:
        trna_frame = pd.read_csv(raw_trnas, sep='\t', skiprows=[0, 2])
        trna_frame.columns = [i.strip() for i in trna_frame.columns]
        # if begin.1 or end.1 are in trnas then drop, else drop the second begin or end
        if 'Begin.1' in trna_frame.columns:
            trna_frame = trna_frame.drop(['Begin.1'], axis=1)
        if 'End.1' in trna_frame.columns:
            trna_frame = trna_frame.drop(['End.1'], axis=1)
        trna_frame = trna_frame.loc[:, get_dups(trna_frame.columns)]
        trna_frame.insert(0, 'fasta', fasta_name)
        return trna_frame
    else:
        logger.warning('No tRNAs were detected, no trnas.tsv file will be created.')
        return None


RAW_RRNA_COLUMNS = ['scaffold', 'tool_name', 'type', 'begin', 'end', 'e-value', 'strand', 'empty', 'note']
RRNA_COLUMNS = ['fasta', 'begin', 'end', 'strand', 'type', 'e-value', 'note']


def run_barrnap(fasta, fasta_name, logger, threads=10, verbose=True):
    raw_rrna_str = run_process(['barrnap', '--threads', str(threads), fasta], logger, capture_stdout=True, check=False,
                               verbose=verbose)
    raw_rrna_table = pd.read_csv(io.StringIO(raw_rrna_str), skiprows=1, sep='\t', header=None,
                                 names=RAW_RRNA_COLUMNS, index_col=0)
    rrna_table_rows = list()
    for gene, row in raw_rrna_table.iterrows():
        rrna_row_dict = {entry.split('=')[0]: entry.split('=')[1] for entry in row['note'].split(';')}
        rrna_table_rows.append([fasta_name, row.begin, row.end, row.strand, rrna_row_dict['Name'].replace('_', ' '),
                                row['e-value'], rrna_row_dict.get('note', '')])
    if len(raw_rrna_table) > 0:
        return pd.DataFrame(rrna_table_rows, index=raw_rrna_table.index, columns=RRNA_COLUMNS).reset_index()
    else:
        logger.warning('No rRNAs were detected, no rrnas.tsv file will be created.')
        return None


def make_trnas_interval(scaffold, row, i):
    if row['Begin'] < row['End']:
        begin = row['Begin']
        end = row['End']
        strand = '+'
    else:
        begin = row['End']
        end = row['Begin']
        strand = '-'
    metadata = {'source': 'tRNAscan-SE', 'type': 'tRNA', 'score': row['Score'], 'strand': strand, 'phase': 0,
                'ID': '%s_tRNA_%s' % (scaffold, i), 'codon': row['Codon'], 'product': 'tRNA-%s' % row['Type']}
    if not pd.isna(row['Note']):
        metadata['Note'] = row['Note']
    return begin, end, metadata


def make_rrnas_interval(scaffold, row, i):
    metadata = {'source': 'barrnap', 'type': 'rRNA', 'score': row['e-value'], 'strand': row['strand'], 'phase': 0,
                'ID': '%s_rRNA_%s' % (scaffold, i), 'product': '%s ribosomal RNA' % row['type'].split()[0],
                'gene': row['type']}
    if not pd.isna(row['note']):
        metadata['Note'] = row['note']
    return row['begin'], row['end'], metadata


# TODO: make it take an input and output gff location and not overwrite
# TODO: for some reason 1 is getting added to intervals when added to gff
def add_intervals_to_gff(annotations_loc, gff_loc, len_dict, interval_function, groupby_column, logger):
    # get fasta length dict so we can merge, I'd love to be able to get this from somewhere else
    annotation_frame = pd.read_csv(annotations_loc, sep='\t')
    # process trnas to intervals
    annotation_dict = dict()
    for scaffold, frame in annotation_frame.groupby(groupby_column):
        if type(scaffold) is not str:
            scaffold = str(scaffold)
        else:
            scaffold = scaffold.strip()
        im = IntervalMetadata(len_dict[scaffold])
        for i, (_, row) in enumerate(frame.iterrows()):
            i += 1
            begin, end, metadata = interval_function(scaffold, row, i)
            begin -= 1
            if begin < 0:
                begin = 0
                logger.warning('Interval added to gff started less than zero, set to zero')
            im.add(bounds=[(begin, end)], metadata=metadata)
        annotation_dict[scaffold] = im
    # add trna intervals to gff
    gff = list(read_sequence(gff_loc, format='gff3'))
    with open(gff_loc, 'w') as f:
        for scaffold, gff_intervals in gff:
            gff_intervals = IntervalMetadata(len_dict[scaffold], gff_intervals)
            if scaffold in annotation_dict:
                gff_intervals.merge(annotation_dict[scaffold])
                gff_intervals.sort()
            f.write(gff_intervals.write(io.StringIO(), format='gff3', seq_id=scaffold).getvalue())


def do_blast_style_search(query_db, target_db, working_dir, db_handler, formater, logger,
                          db_name='database', bit_score_threshold=60, rbh_bit_score_threshold=350, threads=10,
                          verbose=False):
    """A convenience function to do a blast style reciprocal best hits search"""
    # Get kegg hits
    logger.info('Getting forward best hits from %s' % db_name)
    forward_hits = get_best_hits(query_db, target_db, logger, working_dir, 'gene', db_name, bit_score_threshold,
                                 threads, verbose=verbose)
    if stat(forward_hits).st_size == 0:
        return pd.DataFrame()
    logger.info('Getting reverse best hits from %s' % db_name)
    reverse_hits = get_reciprocal_best_hits(query_db, target_db, logger, working_dir, 'gene', db_name,
                                            bit_score_threshold, rbh_bit_score_threshold, threads, verbose=verbose)
    hits = process_reciprocal_best_hits(forward_hits, reverse_hits, db_name)
    logger.info('Getting descriptions of hits from %s' % (db_name))
    if '%s_description' % db_name in db_handler.get_database_names():
        header_dict = db_handler.get_descriptions(hits['%s_hit' % db_name], '%s_description' % db_name)
    else:
        header_dict = multigrep(hits[f'{db_name}_hit'], f'{target_db}_h', logger, '\x00', working_dir)
    hits = formater(hits, header_dict)
    return hits


def count_motifs(gene_faa, motif='(C..CH)'):
    motif_count_dict = dict()
    for seq in read_sequence(gene_faa, format='fasta'):
        motif_count_dict[seq.metadata['id']] = len(list(seq.find_with_regex(motif)))
    return motif_count_dict


def strip_endings(text, suffixes: list):
    for suffix in suffixes:
        if text.endswith(suffix):
            text = text[:len(text) - len(suffix)]
    return text


def process_custom_dbs(custom_fasta_loc, custom_db_name, output_dir, logger, threads=1, verbose=False):
    # if none is passed from argparse then set to tuple of len 0
    mkdir(output_dir)

    if custom_fasta_loc is None:
        custom_fasta_loc = ()
    if custom_db_name is None:
        custom_db_name = ()
    if len(custom_fasta_loc) != len(custom_db_name):
        raise ValueError('Lengths of custom db fasta list and custom db name list must be the same.')
    custom_dbs = {custom_db_name[i]: custom_fasta_loc[i] for i in range(len(custom_db_name))}
    custom_db_locs = dict()
    for db_name, db_loc in custom_dbs.items():
        custom_db_loc = path.join(output_dir, '%s.custom.mmsdb' % db_name)
        make_mmseqs_db(db_loc, custom_db_loc, logger, threads=threads, verbose=verbose)
        custom_db_locs[db_name] = custom_db_loc
    return custom_db_locs


def process_custom_hmms(custom_hmm_loc, custom_hmm_name, logger, verbose=False):
    if custom_hmm_loc is None:
        custom_hmm_loc = ()
    if custom_hmm_name is None:
        custom_hmm_name = ()
    if len(custom_hmm_loc) != len(custom_hmm_name):
        raise ValueError('Lengths of custom db hmm list and custom hmm db name list must be the same.')
    custom_hmm_locs = dict()
    for i in range(len(custom_hmm_name)):
        run_process(['hmmpress', '-f', custom_hmm_loc[i]], logger, verbose=verbose)  # all are pressed just in case
        custom_hmm_locs[custom_hmm_name[i]] = custom_hmm_loc[i]
    return custom_hmm_locs

def process_custom_hmm_cutoffs(custom_hmm_cutoffs_loc, custom_hmm_name, logger:logging.Logger,  verbose=False):
    if custom_hmm_cutoffs_loc is None:
        return {}
    if custom_hmm_name is None:
        raise ValueError("You can't use the custom_hmm_cutoffs_loc argument without the custom_hmm_name and"
                         " custom_hmm_locs aguments specified.")
    if len(custom_hmm_cutoffs_loc) != len(custom_hmm_name):
        logger.warning(f"Custom hmm cutoffs and descriptions were only provided to the first {len(custom_hmm_cutoffs_loc)}."
                      " The rest of the custom hmms will use standard cutoffs and have no descriptions.")
    return {custom_hmm_name[i]:j for i, j in enumerate(custom_hmm_cutoffs_loc)}


class Annotation:
    def __init__(self, name, scaffolds, genes_faa, genes_fna, gff, gbk, annotations, trnas, rrnas):
        # TODO: get abspath for every input file/dir
        # TODO: check that files exist
        self.name = name
        self.scaffolds_loc = scaffolds
        self.genes_faa_loc = genes_faa
        self.genes_fna_loc = genes_fna
        self.gff_loc = gff
        self.gbk_loc = gbk
        self.annotations_loc = annotations
        self.trnas_loc = trnas
        self.rrnas_loc = rrnas

    def get_annotations(self):
        return pd.read_csv(self.annotations_loc, index_col=0, sep='\t')

    def get_trnas(self):
        return pd.read_csv(self.trnas_loc, sep='\t')

    def get_rrnas(self):
        return pd.read_csv(self.rrnas_loc, sep='\t')

def annotate_orfs(gene_faa, db_handler, tmp_dir, logger, custom_db_locs=(), custom_hmm_locs=(),
                  custom_hmm_cutoffs_locs=(), bit_score_threshold=60, rbh_bit_score_threshold=350,
                  kofam_use_dbcan2_thresholds=False, threads=10, verbose=False):
    # run reciprocal best hits searches
    logger.info('Turning genes from prodigal to mmseqs2 db')
    query_db = path.join(tmp_dir, 'gene.mmsdb')
    make_mmseqs_db(gene_faa, query_db, logger, create_index=True, threads=threads, verbose=verbose)

    annotation_list = list()

    if db_handler.config['search_databases'].get('kegg') is not None:
        #TODO Change the get_kegg_description name in function do_blast_style_search to formater
        #TODO think about how this can be consitent with blast and mmseqs
        annotation_list.append(do_blast_style_search(query_db, db_handler.config['search_databases']['kegg'], tmp_dir,
                                                     db_handler, get_kegg_description, logger,
                                                     'kegg', bit_score_threshold, rbh_bit_score_threshold, threads,
                                                     verbose))
    elif db_handler.config['search_databases'].get('kofam_hmm') is not None and db_handler.config['search_databases'].get('kofam_ko_list') is not None:
        logger.info('Getting hits from kofam')
        annotation_list.append(run_hmmscan(genes_faa=gene_faa,
                                           db_loc=db_handler.config['search_databases']['kofam_hmm'],
                                           db_name='kofam_hmm',
                                           output_loc=tmp_dir, #check_impliments
                                           threads=threads, #check_impliments
                                           verbose=verbose,
                                           formater=partial(
                                               kofam_hmmscan_formater,
                                               hmm_info_path=db_handler.config['search_databases']['kofam_ko_list'],
                                               top_hit=True,
                                               use_dbcan2_thresholds=kofam_use_dbcan2_thresholds
                                           ),
                                           logger=logger))
    else:
        logger.warning('No KEGG source provided so distillation will be of limited use.')

    # Get uniref hits
    if db_handler.config['search_databases'].get('uniref') is not None:
        annotation_list.append(do_blast_style_search(query_db, db_handler.config['search_databases']['uniref'], tmp_dir,
                                                     db_handler, get_uniref_description,
                                                     logger, 'uniref', bit_score_threshold,
                                                     rbh_bit_score_threshold, threads, verbose))

    # Get viral hits
    if db_handler.config['search_databases'].get('viral') is not None:
        get_viral_description = partial(get_basic_description, db_name='viral')
        annotation_list.append(do_blast_style_search(query_db, db_handler.config['search_databases']['viral'], tmp_dir,
                                                     db_handler, get_viral_description,
                                                     logger, 'viral', bit_score_threshold,
                                                     rbh_bit_score_threshold, threads, verbose))

    # Get peptidase hits
    if db_handler.config['search_databases'].get('peptidase') is not None:
        annotation_list.append(do_blast_style_search(query_db, db_handler.config['search_databases']['peptidase'], tmp_dir,
                                                     db_handler, get_peptidase_description,
                                                     logger, 'peptidase', bit_score_threshold,
                                                     rbh_bit_score_threshold, threads, verbose))

    # Get pfam hits
    if db_handler.config['search_databases'].get('pfam') is not None:
        logger.info('Getting hits from pfam')
        annotation_list.append(run_mmseqs_profile_search(query_db, db_handler.config['search_databases']['pfam'], tmp_dir,
                                                         logger, output_prefix='pfam', db_handler=db_handler, threads=threads,
                                                         verbose=verbose))

    # use hmmer to detect cazy ids using dbCAN
    if db_handler.config['search_databases'].get('dbcan') is not None:
        logger.info('Getting hits from dbCAN')
        annotation_list.append(run_hmmscan(genes_faa=gene_faa,
                                           db_loc=db_handler.config['search_databases']['dbcan'],
                                           db_name='cazy',
                                           output_loc=tmp_dir,
                                           threads=threads,
                                           formater=partial(
                                               dbcan_hmmscan_formater,
                                               db_name='cazy',
                                               db_handler=db_handler
                                           ),
                                           logger=logger))

    # use hmmer to detect vogdbs
    if db_handler.config['search_databases'].get('vogdb') is not None:
        logger.info('Getting hits from VOGDB')
        annotation_list.append(run_hmmscan(genes_faa=gene_faa,
                                           db_loc=db_handler.config['search_databases']['vogdb'],
                                           db_name='vogdb',
                                           threads=threads,
                                           output_loc=tmp_dir,
                                           formater=partial(
                                               vogdb_hmmscan_formater,
                                               db_name='vogdb',
                                               logger=logger,
                                               db_handler=db_handler
                                           ),
                                           logger=logger))
    for db_name, db_loc in custom_db_locs.items():
        logger.info('Getting hits from %s' % db_name)
        get_custom_description = partial(get_basic_description, db_name=db_name)
        annotation_list.append(do_blast_style_search(query_db, db_loc, tmp_dir, db_handler,
                                                     get_custom_description, logger, db_name,
                                                     bit_score_threshold, rbh_bit_score_threshold, threads,
                                                     verbose))

    # get hits to hmm style custom databases
    for hmm_name, hmm_loc in custom_hmm_locs.items():
        annotation_list.append(run_hmmscan(genes_faa=gene_faa,
                                           db_loc=hmm_loc,
                                           db_name=hmm_name,
                                           threads=threads,
                                           output_loc=tmp_dir,
                                           formater=partial(
                                               generic_hmmscan_formater,
                                               db_name=hmm_name,
                                               hmm_info_path=custom_hmm_cutoffs_locs.get(hmm_name),
                                               top_hit=True
                                           ),
                                           logger=logger))

    annotation_list.append(pd.DataFrame(count_motifs(gene_faa, '(C..CH)'), index=['heme_regulatory_motif_count']).T)

    # merge dataframes
    logger.info('Merging ORF annotations')
    annotations = pd.concat(annotation_list, axis=1, sort=False)

    # get scaffold data and assign grades
    grades = assign_grades(annotations)
    annotations = pd.concat([grades, annotations], axis=1, sort=False)

    return annotations


def annotate_fasta(fasta_loc, fasta_name, output_dir, db_handler, logger, min_contig_size=5000, prodigal_mode='meta',
                   trans_table='11', custom_db_locs=(), custom_hmm_locs=(), custom_hmm_cutoffs_locs=(),
                   bit_score_threshold=60, rbh_bit_score_threshold=350, kofam_use_dbcan2_thresholds=False,
                   skip_trnascan=False, threads=1, rename_bins=True, keep_tmp_dir=False,
                   verbose=False):
    """Annotated a single multifasta file, all file based outputs will be in output_dir"""
    # make temporary directory
    tmp_dir = path.join(output_dir, 'tmp')
    mkdir(tmp_dir)

    # filter input fasta
    filtered_fasta = path.join(tmp_dir, 'filtered_fasta.fa')
    filter_fasta(fasta_loc, min_contig_size, filtered_fasta)

    if stat(filtered_fasta).st_size == 0:
        logger.warning('No sequences were longer than min_contig_size')
        return None
    # predict ORFs with prodigal
    # TODO: handle when prodigal returns no genes
    gene_gff, gene_fna, gene_faa = run_prodigal(filtered_fasta, tmp_dir, logger, mode=prodigal_mode,
                                                trans_table=trans_table, verbose=verbose)

    # annotate ORFs
    annotations = annotate_orfs(gene_faa, db_handler, tmp_dir, logger, custom_db_locs, custom_hmm_locs,
                                custom_hmm_cutoffs_locs, bit_score_threshold, rbh_bit_score_threshold,
                                kofam_use_dbcan2_thresholds, threads, verbose)
    annotations = pd.concat([get_gene_data(gene_faa), annotations], axis=1, sort=False)

    renamed_scaffolds = path.join(output_dir, 'scaffolds.annotated.fa')
    if rename_bins:
        # rename scaffolds to match prodigal names
        prefix = fasta_name
        rename_fasta(filtered_fasta, renamed_scaffolds, prefix=prefix)
    else:
        prefix = None
        copy2(filtered_fasta, renamed_scaffolds)

    # generate fna and faa output files with annotations
    annotated_fna = path.join(output_dir, 'genes.annotated.fna')
    create_annotated_fasta(gene_fna, annotations, annotated_fna, name=prefix)
    annotated_faa = path.join(output_dir, 'genes.annotated.faa')
    create_annotated_fasta(gene_faa, annotations, annotated_faa, name=prefix)

    # rename gff entries to match prodigal names
    renamed_gffs = path.join(output_dir, 'genes.annotated.gff')
    annotate_gff(gene_gff, renamed_gffs, annotations, prefix=prefix)

    # add fasta name to frame and index, append to list
    annotations.insert(0, 'fasta', fasta_name)
    if rename_bins:
        annotations.index = annotations.fasta + '_' + annotations.index
    annotations_loc = path.join(output_dir, 'annotations.tsv')
    annotations.to_csv(annotations_loc, sep='\t')

    # get tRNAs and rRNAs
    len_dict = {i.metadata['id']: len(i) for i in read_sequence(renamed_scaffolds, format='fasta')}
    if not skip_trnascan:
        trna_table = run_trna_scan(renamed_scaffolds, tmp_dir, fasta_name, logger, threads=threads, verbose=verbose)
        if trna_table is not None:
            trna_loc = path.join(output_dir, 'trnas.tsv')
            trna_table.to_csv(trna_loc, sep='\t', index=False)
            add_intervals_to_gff(trna_loc, renamed_gffs, len_dict, make_trnas_interval, 'Name', logger)
        else:
            trna_loc = None
    else:
        trna_loc = None

    rrna_table = run_barrnap(renamed_scaffolds, fasta_name, logger,  threads=threads, verbose=verbose)
    if rrna_table is not None:
        rrna_loc = path.join(output_dir, 'rrnas.tsv')
        rrna_table.to_csv(rrna_loc, sep='\t', index=False)
        add_intervals_to_gff(rrna_loc, renamed_gffs, len_dict, make_rrnas_interval, 'scaffold', logger)
    else:
        rrna_loc = None

    # make genbank file
    current_gbk = path.join(output_dir, '%s.gbk' % fasta_name)
    make_gbk_from_gff_and_fasta(renamed_gffs, renamed_scaffolds, annotated_faa, current_gbk)

    if not keep_tmp_dir:
        rmtree(tmp_dir)

    return Annotation(name=fasta_name, scaffolds=renamed_scaffolds, genes_faa=annotated_faa, genes_fna=annotated_fna,
                      gff=renamed_gffs, gbk=current_gbk, annotations=annotations_loc, trnas=trna_loc, rrnas=rrna_loc)


def get_fasta_name(fasta_loc):
    return path.splitext(path.basename(remove_suffix(fasta_loc, '.gz')))[0]


def annotate_fastas(fasta_locs, output_dir, db_handler, logger, min_contig_size=5000, prodigal_mode='meta',
                    trans_table='11', bit_score_threshold=60, rbh_bit_score_threshold=350, custom_db_name=(),
                    custom_fasta_loc=(), custom_hmm_name=(), custom_hmm_loc=(), custom_hmm_cutoffs_loc=(),
                    kofam_use_dbcan2_thresholds=False, skip_trnascan=False, rename_bins=True, keep_tmp_dir=True,
                    threads=10, verbose=True):
    # check for no conflicting options/configurations
    tmp_dir = path.join(output_dir, 'working_dir')
    mkdir(tmp_dir)

    # setup custom databases to be searched
    custom_db_locs = process_custom_dbs(custom_fasta_loc, custom_db_name, path.join(tmp_dir, 'custom_dbs'),
                                        logger, threads, verbose)
    custom_hmm_locs = process_custom_hmms(custom_hmm_loc, custom_hmm_name, logger)
    custom_hmm_cutoffs_locs= process_custom_hmm_cutoffs(custom_hmm_cutoffs_loc, custom_hmm_name, logger)
    logger.info('Retrieved database locations and descriptions')

    # iterate over list of fastas and annotate each individually
    annotations_list = list()
    for fasta_loc in fasta_locs:
        # get name of file e.g. /home/shaffemi/my_genome.fa -> my_genome
        fasta_name = get_fasta_name(fasta_loc)
        logger.info('Annotating %s' % fasta_name)
        fasta_dir = path.join(tmp_dir, fasta_name)
        mkdir(fasta_dir)
        annotations_list.append(
            annotate_fasta(fasta_loc, fasta_name, fasta_dir, db_handler, logger, min_contig_size, prodigal_mode,
                           trans_table, custom_db_locs, custom_hmm_locs, custom_hmm_cutoffs_locs,
                           bit_score_threshold, rbh_bit_score_threshold, kofam_use_dbcan2_thresholds,
                           skip_trnascan, threads, rename_bins, keep_tmp_dir, verbose))
    logger.info('Annotations complete, processing annotations')

    all_annotations = merge_annotations(annotations_list, output_dir)

    # clean up
    if not keep_tmp_dir:
        rmtree(tmp_dir)
    return all_annotations


def make_fasta_namses_df(fasta_loc):
    fasta_name = get_fasta_name(fasta_loc)
    names = [{'fasta': fasta_name, 'seq': seq.metadata['id']} for seq in read_sequence(fasta_loc,format='fasta')]
    return pd.DataFrame(names)


# TODO: Add force flag to remove output dir if it already exists
# TODO: Add continute flag to continue if output directory already exists
# TODO: make fasta loc either a string or list to remove annotate_bins_cmd and annotate_called_genes_cmd?
def annotate_bins(input_fasta:list, output_dir='.', min_contig_size=2500, prodigal_mode='meta', trans_table='11',
                  bit_score_threshold=60, rbh_bit_score_threshold=350, custom_db_name=(), custom_fasta_loc=(),
                  custom_hmm_name=(), custom_hmm_loc=(), custom_hmm_cutoffs_loc=(), use_uniref=False,
                  use_vogdb=False, kofam_use_dbcan2_thresholds=False,
                  skip_trnascan=False, gtdb_taxonomy=(), checkm_quality=(),
                  rename_bins=True, keep_tmp_dir=True, low_mem_mode=False, threads=10, verbose=True,
                  log_file_path:str=None, join_fastas:bool=False, config_loc:str=None):
    rename_bins = True
    fasta_locs = [j for i in input_fasta for j in glob(i)]
    mkdir(output_dir)
    if log_file_path is None:
        log_file_path = path.join(output_dir, "annotate.log")
    logger = logging.getLogger('annotation_log')
    setup_logger(logger, log_file_path)
    logger.info(f"The log file is created at {log_file_path}.")

    if len(fasta_locs) == 0:
        raise ValueError('Given fasta locations return no paths: %s' % input_fasta)
    fasta_names = [get_fasta_name(i) for i in fasta_locs]
    if len(fasta_names) != len(set(fasta_names)):
        raise ValueError('Genome file names must be unique. At least one name appears twice in this search.')
    logger.info('%s FASTAs found' % len(fasta_locs))
    # set up
    db_handler = DatabaseHandler(logger, config_loc)
    db_handler.filter_db_locs(low_mem_mode, use_uniref,
                              use_vogdb, master_list=MAG_DBS_TO_ANNOTATE,)
    db_conf = db_handler.get_settings_str()
    logger.info(f"Starting the Annotation of Bins with database configuration: \n {db_conf}")

    # check inputs
    prodigal_modes = ['train', 'meta', 'single']
    if prodigal_mode not in prodigal_modes:
        raise ValueError('Prodigal mode must be one of %s.' % ', '.join(prodigal_modes))
    elif prodigal_mode in ['normal', 'single']:
        logger.warning('When running prodigal in single mode your bins must have long contigs (average length >3 Kbp), '
                      'be long enough (total length > 500 Kbp) and have very low contamination in order for prodigal '
                      'training to work well.')

    # prodigal_trans_tables = ['auto'] + [str(i) for i in range(1, 26)]
    prodigal_trans_tables = [str(i) for i in range(1, 26)]
    if trans_table not in prodigal_trans_tables:
        # raise ValueError('Prodigal translation table must be 1-25 or auto')
        raise ValueError('Prodigal translation table must be 1-25')


    all_annotations = annotate_fastas(fasta_locs, output_dir, db_handler, logger, min_contig_size, prodigal_mode, trans_table,
                                      bit_score_threshold, rbh_bit_score_threshold, custom_db_name, custom_fasta_loc,
                                      custom_hmm_name, custom_hmm_loc, custom_hmm_cutoffs_loc,
                                      kofam_use_dbcan2_thresholds, skip_trnascan, rename_bins, keep_tmp_dir,
                                      threads, verbose)
    # if given add taxonomy information
    if len(gtdb_taxonomy) > 0:
        gtdb_taxonomy = pd.concat([pd.read_csv(i, sep='\t', index_col=0) for i in gtdb_taxonomy])
        taxonomy = list()
        taxonomy_missing_bins = list()
        for i in all_annotations.fasta:
            # add taxonomy
            if i in gtdb_taxonomy.index:
                taxonomy.append(gtdb_taxonomy.loc[i, 'classification'])
            else:
                taxonomy.append(i)
                taxonomy_missing_bins.append(i)
        for i in set(taxonomy_missing_bins):
            logger.warning('Bin %s was not found in taxonomy file, replaced with bin name.' % i)
        all_annotations['bin_taxonomy'] = taxonomy
    # if given add quality information
    if len(checkm_quality) > 0:
        checkm_quality = pd.concat([pd.read_csv(i, sep='\t', index_col=0) for i in checkm_quality])
        checkm_quality.index = [strip_endings(i, ['.fa', '.fasta', '.fna']) for i in checkm_quality.index]

        completeness = list()
        contamination = list()
        quality_missing_bins = list()
        for i in all_annotations.fasta:
            # add completeness and contamination
            if i in checkm_quality.index:
                completeness.append(checkm_quality.loc[i, 'Completeness'])
                contamination.append(checkm_quality.loc[i, 'Contamination'])
            else:
                completeness.append(0)
                contamination.append(100)
                quality_missing_bins.append(i)
            for j in set(quality_missing_bins):
                logger.warning('Bin %s was not found in quality file, '
                              'replaced with completeness 0 and contamination 100.' % j)
        all_annotations['bin_completeness'] = completeness
        all_annotations['bin_contamination'] = contamination
    all_annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')

    logger.info("Completed annotations")


def annotate_called_genes_cmd(input_faa, output_dir='.', bit_score_threshold=60,
                              rbh_bit_score_threshold=350,
                              custom_db_name=(), custom_fasta_loc=(), custom_hmm_loc=(), custom_hmm_name=(),
                              custom_hmm_cutoffs_loc=(), use_uniref=False, log_file_path:str=None,
                              use_vogdb=False, kofam_use_dbcan2_thresholds=False, rename_genes=True,
                              keep_tmp_dir=True, low_mem_mode=False, threads=10, verbose=True,
                              config_loc:str=None):
    fasta_locs = glob(input_faa)
    annotate_called_genes(fasta_locs, output_dir, bit_score_threshold, rbh_bit_score_threshold,
                          custom_db_name, custom_fasta_loc, custom_hmm_loc, custom_hmm_name,
                          custom_hmm_cutoffs_loc, use_uniref, use_vogdb, kofam_use_dbcan2_thresholds,
                          rename_genes, keep_tmp_dir, low_mem_mode, threads, verbose, None,
                          config_loc)

def perform_fasta_checks(fasta_locs, logger):
    """Perform all checks related to integrity of fastas"""
    if len(fasta_locs) == 0:
        raise ValueError('Given fasta locations returns no paths: %s' % input_faa)
    logger.info('%s fastas found' % len(fasta_locs))
    logger.info("Checking for duplicate names")
    names = [get_fasta_name(i) for i in fasta_locs]
    try:
        if len(names) != len(set(names)):
            raise ValueError('Genome file names must be unique. At least one name appears twice in this search.')
        fastas_dup_check(fasta_locs, '>')
    except ValueError as error:
        logger.critical(error)
        raise error
    logger.info("No duplicate names found")

def annotate_called_genes(fasta_locs, output_dir='.', bit_score_threshold=60, rbh_bit_score_threshold=350,
                          custom_db_name=(), custom_fasta_loc=(), custom_hmm_loc=(), custom_hmm_name=(),
                          custom_hmm_cutoffs_loc=(), use_uniref=False,
                          use_vogdb=False, kofam_use_dbcan2_thresholds=False, rename_genes=True, keep_tmp_dir=True,
                          low_mem_mode=False, threads=10, verbose=True, log_file_path:str=None, config_loc:str=None):
    mkdir(output_dir)

    # Get a logger
    if log_file_path is None:
        log_file_path = path.join(output_dir, "annotate.log")
    logger = logging.getLogger('annotation_log')
    setup_logger(logger, log_file_path)
    logger.info(f"The log file is created at {log_file_path}")
    # get database locations
    db_handler = DatabaseHandler(logger, config_loc)
    db_handler.filter_db_locs(low_mem_mode, use_uniref, use_vogdb, master_list=MAG_DBS_TO_ANNOTATE)

    # Check fastas
    perform_fasta_checks(fasta_locs, logger)

    logger.info(f"Starting the Annotation of Genes with database configuration: \n {db_handler.get_settings_str()}")

    logger.info("Checking for duplicate names")
    fasta_names = [get_fasta_name(i) for i in fasta_locs]
    if len(fasta_names) != len(set(fasta_names)):
        raise ValueError('Genome file names must be unique. At least one name appears twice in this search.')
    logger.info("No duplicate names found")


    tmp_dir = path.join(output_dir, 'working_dir')
    mkdir(tmp_dir)

    # setup custom databases to be searched
    custom_db_locs = process_custom_dbs(custom_fasta_loc, custom_db_name, path.join(tmp_dir, 'custom_dbs'),
                                        logger, threads, verbose)
    custom_hmm_locs = process_custom_hmms(custom_hmm_loc, custom_hmm_name, logger)
    custom_hmm_cutoffs_locs= process_custom_hmm_cutoffs(custom_hmm_cutoffs_loc, custom_hmm_name, logger)
    logger.info('Retrieved database locations and descriptions')

    # annotate
    annotation_locs = list()
    faa_locs = list()
    for fasta_loc in fasta_locs:
        # set up
        fasta_name = get_fasta_name(fasta_loc)
        fasta_dir = path.join(tmp_dir, fasta_name)
        mkdir(fasta_dir)

        # annotate
        annotations = annotate_orfs(fasta_loc, db_handler, fasta_dir, logger, custom_db_locs, custom_hmm_locs,
                                    custom_hmm_cutoffs_locs, bit_score_threshold, rbh_bit_score_threshold,
                                    kofam_use_dbcan2_thresholds, threads, verbose)

        annotated_faa = path.join(fasta_dir, 'genes.faa')

        create_annotated_fasta(fasta_loc, annotations, annotated_faa, name=fasta_name)
        faa_locs.append(annotated_faa)

        # add fasta name to frame and index, write file
        annotations.insert(0, 'fasta', fasta_name)
        if rename_genes:
            annotations.index = annotations.fasta + '_' + annotations.index
        annotation_loc = path.join(fasta_dir, 'annotations.tsv')
        annotations.to_csv(annotation_loc, sep='\t')
        annotation_locs.append(annotation_loc)

    # merge
    all_annotations = pd.concat([pd.read_csv(i, sep='\t', index_col=0) for i in annotation_locs], sort=False)
    all_annotations = all_annotations.sort_values('fasta')
    all_annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')
    merge_files(faa_locs, path.join(output_dir, 'genes.faa'))

    # clean up
    if not keep_tmp_dir:
        rmtree(tmp_dir)

    logger.info("Completed annotations")


def merge_annotations(annotations_list, output_dir, write_annotations=False):

    # merge annotation dicts
    all_annotations = pd.concat([i.get_annotations() for i in annotations_list if i is not None], sort=False)
    all_annotations = all_annotations.sort_values(['fasta', 'scaffold', 'gene_position'])

    # merge gene files
    merge_files([i.genes_fna_loc for i in annotations_list if i is not None], path.join(output_dir, 'genes.fna'))
    merge_files([i.genes_faa_loc for i in annotations_list if i is not None], path.join(output_dir, 'genes.faa'))
    merge_files([i.scaffolds_loc for i in annotations_list if i is not None], path.join(output_dir, 'scaffolds.fna'))
    merge_files([i.gff_loc for i in annotations_list if i is not None], path.join(output_dir, 'genes.gff'), True)
    trnas_locs = [i.trnas_loc for i in annotations_list if i is not None if i.trnas_loc is not None]
    if len(trnas_locs) > 0:
        merge_files(trnas_locs, path.join(output_dir, 'trnas.tsv'), True)
    rrnas_locs = [i.rrnas_loc for i in annotations_list if i is not None if i.rrnas_loc is not None]
    if len(rrnas_locs) > 0:
        merge_files(rrnas_locs, path.join(output_dir, 'rrnas.tsv'), True)

    # make output gbk dir
    gbk_dir = path.join(output_dir, 'genbank')
    mkdir(gbk_dir)
    for anno in annotations_list:
        if anno is not None:
            # TODO: make annotate_fasta generate a genbank dir and then copy it's contents, get rid of Annotation.name
            if path.isfile(anno.gbk_loc):
                copy2(anno.gbk_loc, path.join(gbk_dir, '%s.gbk' % anno.name))
            else:
                for gbk_loc in glob(path.join(anno.gbk_loc, '*.gbk')):
                    copy2(gbk_loc, gbk_dir)
    if write_annotations:
        all_annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')
    else:
        return all_annotations


def merge_annotations_cmd(input_dirs, output_dir):

    mkdir(output_dir)
    # Get a logger
    annotations_list = list()
    log_file_path = path.join(output_dir, "annotate.log")
    logger = logging.getLogger('annotate.log')
    setup_logger(logger, log_file_path)
    logger.info(f"The log file is created at {log_file_path}")

    # make Annotation objects per directory
    for i, annotation_dir in enumerate(glob(input_dirs)):
        if not path.exists(path.join(annotation_dir, "annotations.tsv")):
            logger.warning("Skipping the fallowing directory because"
                           f" no 'annotations.tsv' file was found: {annotation_dir}" )
            continue
        name = 'annotation_%s' % i
        scaffolds = path.join(annotation_dir, 'scaffolds.fna')
        genes_faa = path.join(annotation_dir, 'genes.faa')
        genes_fna = path.join(annotation_dir, 'genes.fna')
        gff = path.join(annotation_dir, 'genes.gff')
        gbk = path.join(annotation_dir, 'genbank')
        annotations = path.join(annotation_dir, 'annotations.tsv')
        trnas = path.join(annotation_dir, 'trnas.tsv')
        if not path.isfile(trnas):
            logger.warning("No trnas.tsv file found in directory %s" % annotation_dir)
            trnas = None
        rrnas = path.join(annotation_dir, 'rrnas.tsv')
        if not path.isfile(rrnas):
            logger.warning("No rrnas.tsv file found in directory %s" % annotation_dir)
            rrnas = None
        annotations_list.append(Annotation(name=name, scaffolds=scaffolds, genes_faa=genes_faa, genes_fna=genes_fna,
                                           gff=gff, gbk=gbk, annotations=annotations, trnas=trnas, rrnas=rrnas))
    # run merge_annotations
    merge_annotations(annotations_list, output_dir, write_annotations=True)

