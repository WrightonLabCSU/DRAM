import pandas as pd
from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.utils import get_database_locs
from mag_annotator.summarize_vgfs import get_ids_from_row, filter_to_amgs

# TODO: filter by taxonomic level, completeness, contamination
# TODO: filter scaffolds file, gff file
# TODO: add negate, aka pull not from given list


def pull_sequences(input_annotations, input_fasta, output_fasta, fastas=None, scaffolds=None, genes=None,
                   identifiers=None, categories=None, taxonomy=None, completeness=None, contamination=None,
                   amg_flags=None, aux_scores=None, virsorter_category=None, putative_amgs=False, max_auxiliary_score=3,
                   remove_transposons=False, remove_fs=False, remove_js=False):
    annotations = pd.read_csv(input_annotations, sep='\t', index_col=0)

    # first filter based on specific names
    specific_genes_to_keep = list()
    # filter fastas
    if fastas is not None:
        for fasta in fastas:
            specific_genes_to_keep += list(annotations.loc[annotations['fasta'] == fasta].index)
    # filter scaffolds
    if scaffolds is not None:
        for scaffold in scaffolds:
            specific_genes_to_keep += list(annotations.loc[annotations['scaffold'] == scaffold].index)
    # filter genes
    if genes is not None:
        specific_genes_to_keep += genes
    # filter down annotations based on specific genes
    if len(specific_genes_to_keep) > 0:
        annotations = annotations.loc[specific_genes_to_keep]

    # filter based on annotations
    if (identifiers is not None) or (categories is not None):
        annotation_genes_to_keep = list()
        gene_to_ids = dict()
        for i, row in annotations.iterrows():
            row_ids = get_ids_from_row(row)
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
            db_locs = get_database_locs()
            genome_summary_form = pd.read_csv(db_locs['genome_summary_form'], sep='\t')
            for level in ['module', 'sheet', 'header', 'subheader']:
                for category, frame in genome_summary_form.loc[~pd.isna(genome_summary_form[level])].groupby(level):
                    if category in categories:
                        for gene, ids in gene_to_ids.items():
                            if len(ids & set(frame['gene_id'])) > 0:
                                annotation_genes_to_keep.append(gene)
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
                                     remove_transposons=remove_transposons, remove_fs=remove_fs, remove_js=remove_js)
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
