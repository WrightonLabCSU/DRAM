import pandas as pd
from skbio import read as read_sequence
from skbio import write as write_sequence
from os import mkdir, path
import warnings

from mag_annotator.utils import get_genes_from_identifiers
from mag_annotator.summarize_vgfs import filter_to_amgs

# TODO: filter by taxonomic level, completeness, contamination
# TODO: filter scaffolds file, gff file
# TODO: add negate, aka pull not from given list


# def pull_sequences(input_annotations, input_fasta, output_fasta, fastas=None, scaffolds=None, genes=None,
#                    identifiers=None, categories=None, taxonomy=None, completeness=None, contamination=None,
#                    amg_flags=None, aux_scores=None, virsorter_category=None, putative_amgs=False, max_auxiliary_score=3,
#                    remove_transposons=False, remove_fs=False, remove_js=False):
def pull_sequences(input_annotations, input_fasta, output_fasta, fastas=None, scaffolds=None, genes=None,
                   identifiers=None, categories=None, taxonomy=None, completeness=None, contamination=None,
                   amg_flags=None, aux_scores=None, virsorter_category=None, putative_amgs=False,
                   max_auxiliary_score=3, remove_transposons=False, remove_fs=False):
    annotations = pd.read_csv(input_annotations, sep='\t', index_col=0)

    # first filter based on specific names
    annotation_genes_to_keep = get_genes_from_identifiers(annotations, genes, fastas, scaffolds, identifiers,
                                                          categories)
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
        # annotations = filter_to_amgs(annotations.fillna(''), max_aux=max_auxiliary_score,
        #                              remove_transposons=remove_transposons, remove_fs=remove_fs, remove_js=remove_js)
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
                           scaffolds_loc=None, distance_genes=None, distance_bp=None):
    # check inputs, make output
    if distance_genes is None and distance_bp is None:
        raise ValueError("Must provide distance away in bp, genes or both.")

    # get data
    annotations = pd.read_csv(input_file, sep='\t', index_col=0)
    genes_from_ids = get_genes_from_identifiers(annotations, genes=genes, identifiers=identifiers,
                                                categories=categories)
    if len(genes_from_ids) == 0:
        raise ValueError("No genes were found based on your filtering parameters. No neighborhoods will be generated.")

    mkdir(output_dir)

    neighborhood_all_annotations = find_neighborhoods(annotations, genes_from_ids, distance_bp, distance_genes)
    neighborhood_all_annotations.to_csv(path.join(output_dir, 'neighborhood_annotations.tsv'), sep='\t')

    # filter files if given
    if genes_loc is not None:
        output_fasta_generator = (i for i in read_sequence(genes_loc, format='fasta')
                                  if i.metadata['id'] in neighborhood_all_annotations.index)
        # TODO: potentially generate one fasta file per neighborhood
        write_sequence(output_fasta_generator, format='fasta',
                       into=path.join(output_dir, 'neighborhood_genes.%s' % genes_loc.split('.')[-1]))

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
