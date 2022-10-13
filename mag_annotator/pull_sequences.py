import pandas as pd
from skbio import read as read_sequence
from skbio import write as write_sequence
from os import mkdir, path
import warnings
import logging

from mag_annotator.summarize_vgfs import filter_to_amgs
from mag_annotator.utils import setup_logger
from mag_annotator.summarize_genomes import get_ids_from_annotations_by_row
from mag_annotator.database_handler import DatabaseHandler

# TODO: filter by taxonomic level, completeness, contamination
# TODO: filter scaffolds file, gff file
# TODO: add negate, aka pull not from given list


def get_genes_from_identifiers(annotations, genes=None, fastas=None, scaffolds=None, identifiers=None, categories=None,
                               custom_distillate=None):
    specific_genes_to_keep = list()
    # filter fastas
    if fastas is not None:
        for fasta in fastas:
            specific_genes_to_keep += list(annotations.loc[annotations['fasta'] == fasta].index)
    # filter scaffolds
    # TODO: remove this functionality or modify since scaffolds are guaranteed unique
    if scaffolds is not None:
        for scaffold in scaffolds:
            specific_genes_to_keep += list(annotations.loc[annotations['scaffold'] == scaffold].index)
    # filter genes
    if genes is not None:
        specific_genes_to_keep += genes
    # filter down annotations based on specific genes
    if len(specific_genes_to_keep) > 0:
        annotations_to_keep = annotations.loc[specific_genes_to_keep]
    else:
        annotations_to_keep = annotations

    # filter based on annotations
    if (identifiers is not None) or (categories is not None):
        annotation_genes_to_keep = list()

        # make a dictionary of genes to
        gene_to_ids = dict()
        for i, row_ids in get_ids_from_annotations_by_row(annotations_to_keep).iteritems():
           if len(row_ids) > 0:
                gene_to_ids[i] = set(row_ids)

        # get genes with ids
        if identifiers is not None:
            identifiers = set(identifiers)
            for gene, ids in gene_to_ids.items():
                if len(set(ids) & set(identifiers)) > 0:
                    annotation_genes_to_keep.append(gene)

        # get genes from distillate categories
        if categories is not None:
            database_handler = DatabaseHandler()
            genome_summary_form = pd.read_csv(database_handler.dram_sheet_locs['genome_summary_form'], sep='\t')
            if custom_distillate is not None:
                genome_summary_form = pd.concat([genome_summary_form, pd.read_csv(custom_distillate, sep='\t')])
            for level in ['module', 'sheet', 'header', 'subheader']:
                for category, frame in genome_summary_form.loc[~pd.isna(genome_summary_form[level])].groupby(level):
                    if category in categories:
                        for gene, ids in gene_to_ids.items():
                            if len(ids & set(frame['gene_id'])) > 0:
                                annotation_genes_to_keep.append(gene)
    else:
        annotation_genes_to_keep = list(annotations_to_keep.index)
    return annotation_genes_to_keep


def pull_sequences(input_tsv, input_fasta, output_fasta, fastas=None, scaffolds=None, genes=None,
                   identifiers=None, adjective_sheet:str=None, categories=None, taxonomy=None, completeness=None, contamination=None,
                   amg_flags=None, aux_scores=None, virsorter_category=None, putative_amgs=False,
                   max_auxiliary_score=3, remove_transposons=False, remove_fs=False, custom_distillate=None):

    #setup logger
    logger = logging.getLogger('strainer_log')
    setup_logger(logger)

    try:
        annotations = pd.read_csv(input_tsv, sep='\t', index_col=0)
        if adjective_sheet is not None:
            adjective_sheet = pd.read_csv(input_tsv, sep='\t', index_col=0) 
            annotations = annotations.loc[adjective_sheet.index]
        annotation_genes_to_keep = get_genes_from_identifiers(annotations, genes, fastas, scaffolds, identifiers,
                                                              categories, custom_distillate)
        annotations = annotations.loc[annotation_genes_to_keep]
        if len(annotations) == 0:
            raise ValueError("Categories or identifiers provided yielded no annotations")


        # DRAM specific filtering
        if taxonomy is not None:
            taxonomy = set(taxonomy)
            annotations = annotations.loc[[len(set(i.split(';')) & taxonomy) > 0 for i in annotations['bin_taxonomy']]]
        if completeness is not None:
            annotations = annotations.loc[annotations['bin_completeness'].astype(float) > completeness]
        if contamination is not None:
            annotations = annotations.loc[annotations['bin_contamination'].astype(float) < contamination]
        if len(annotations) == 0:
            raise ValueError("DRAM filters yielded no annotations")
        # DRAM-v specific filtering
        if putative_amgs:  # get potential amgs
            annotations = filter_to_amgs(annotations.fillna(''), max_aux=max_auxiliary_score,
                                         remove_transposons=remove_transposons, remove_fs=remove_fs)
        else:
            # filter based on virsorter categories
            if virsorter_category is not None:
                annotations = annotations.loc[[i in virsorter_category for i in annotations.virsorter]]
            # filter based on aux scores
            if aux_scores is not None:
                annotations = annotations.loc[[i in aux_scores for i in annotations.auxiliary_score]]
            # filter based on amg flags
            if amg_flags is not None:
                amg_flags = set(amg_flags)
                annotations = annotations.loc[[len(set(i) & amg_flags) > 0 if not pd.isna(i) else False
                                               for i in annotations.amg_flags]]
        if len(annotations) == 0:
            raise ValueError("DRAM-v filters yielded no annotations")

        # make output
        output_fasta_generator = (i for i in read_sequence(input_fasta, format='fasta')
                                  if i.metadata['id'] in annotations.index)
        write_sequence(output_fasta_generator, format='fasta', into=output_fasta)
    except KeyError as err:
        logger.error(err)
        logger.critical("It looks like your input files are not appropriate for the filter you are trying to apply."
                     " For example using a Viral or AMG filter on an annotations that has no viral information.")
        
        raise err
    except Exception as err:
        logger.error(err)
        logger.critical('Unknown error, exiting')
        raise err


def find_neighborhoods(annotations, genes_from_ids, distance_bp=None, distance_genes=None):
    # get neighborhoods as dataframes
    neighborhood_frames = list()
    for neighborhood_number, gene in enumerate(genes_from_ids):
        gene_row = annotations.loc[gene]
        scaffold_annotations = annotations.loc[annotations['scaffold'] == gene_row['scaffold']]
        # get neighbors based on bp
        if distance_bp is not None:
            right_dist = gene_row['end_position'] + distance_bp
            left_dist = gene_row['start_position'] - distance_bp
            neighborhood_annotations = scaffold_annotations.loc[(scaffold_annotations['end_position'] >= left_dist) &
                                                                (scaffold_annotations['start_position'] <= right_dist)]
        else:
            neighborhood_annotations = scaffold_annotations
        # get neighbors based on annotations
        if distance_genes is not None:
            right_genes = gene_row['gene_position'] + distance_genes
            left_genes = gene_row['gene_position'] - distance_genes
            neighborhood_annotations = scaffold_annotations.loc[(scaffold_annotations['gene_position'] >= left_genes) &
                                                                (scaffold_annotations['gene_position'] <= right_genes)]
        # add neighborhood number and neighborhood center as columns
        neighborhood_annotations['neighborhood_number'] = neighborhood_number
        neighborhood_annotations['neighborhood_center'] = [i == gene for i in neighborhood_annotations.index]
        neighborhood_frames.append(neighborhood_annotations)
        if len(neighborhood_annotations) == 0:
            warnings.warn("")

    # merge data frames and write to file
    return pd.concat(neighborhood_frames)


def get_gene_neighborhoods(input_file, output_dir, genes=None, identifiers=None, categories=None, genes_loc=None,
                           scaffolds_loc=None, distance_genes=None, distance_bp=None, custom_distillate=None):
    # check inputs, make output
    if distance_genes is None and distance_bp is None:
        raise ValueError("Must provide distance away in bp, genes or both.")

    # get data
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    genes_from_ids = get_genes_from_identifiers(annotations, genes=genes, identifiers=identifiers,
                                                categories=categories, custom_distillate=custom_distillate)
    if len(genes_from_ids) == 0:
        raise ValueError("No genes were found based on your filtering parameters. No neighborhoods will be generated.")

    mkdir(output_dir)

    #setup logger
    # Get a logger
    # if log_file_path is None:
    log_file_path = path.join(output_dir, "strainer.log")
    logger = logging.getLogger('strainer_log')
    setup_logger(logger, log_file_path)
    logger.info(f"The log file is created at {log_file_path}")

    neighborhood_all_annotations = find_neighborhoods(annotations, genes_from_ids, distance_bp, distance_genes)
    neighborhood_all_annotations.to_csv(path.join(output_dir, 'neighborhood_annotations.tsv'), sep='\t')
    logger.info("Neighborhood Annotations witten to tsv")
    # filter files if given
    if genes_loc is not None:
        output_fasta_generator = (i for i in read_sequence(genes_loc, format='fasta')
                                  if i.metadata['id'] in neighborhood_all_annotations.index)
        # TODO: potentially generate one fasta file per neighborhood
        write_sequence(output_fasta_generator, format='fasta',
                       into=path.join(output_dir, 'neighborhood_genes.%s' % genes_loc.split('.')[-1]))
        logger.info("Gene Neighborhood fasta generated")
    if scaffolds_loc is not None:
        neighborhood_all_annotations['scaffold_mod'] = ['%s_%s' % (row['fasta'], row['scaffold'])
                                                        for i, row in neighborhood_all_annotations.iterrows()]
        neighborhood_scaffolds = list()
        for scaffold in read_sequence(scaffolds_loc, format='fasta'):
            if scaffold.metadata['id'] in neighborhood_all_annotations['scaffold_mod'].values:
                scaffold_frame = neighborhood_all_annotations.loc[neighborhood_all_annotations['scaffold_mod'] ==
                                                                  scaffold.metadata['id']]
                for neighborhood, neighborhood_frame in scaffold_frame.groupby('neighborhood_number'):
                    neighborhood_frame = neighborhood_frame.sort_values('start_position')
                    neighborhood_scaffolds.append(scaffold[neighborhood_frame['start_position'][0]:
                                                           neighborhood_frame['end_position'][-1]])
        write_sequence((i for i in neighborhood_scaffolds), format='fasta',
                       into=path.join(output_dir, 'neighborhood_scaffolds.fna'))
        logger.info("Scaffolds Neighborhood fasta generated")
