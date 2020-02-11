import pandas as pd
from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.utils import get_database_locs
from mag_annotator.summarize_vgfs import get_ids_from_row

# TODO: filter by taxonomic level, completeness, contamination
# TODO: filter scaffolds file, gff file
# TODO: add negate, aka pull not from given list
# TODO: update get_ids_from_row to get PFAM ids
# TODO: pull sequences based on distillate categories

def pull_sequences(input_annotations, input_fasta, output_fasta, fastas=None, scaffolds=None, genes=None,
                   identifiers=None, categories=None):
    annotations = pd.read_csv(input_annotations, sep='\t', index_col=0)

    genes_to_keep = list()

    # filter fastas
    if fastas is not None:
        for fasta in fastas:
            genes_to_keep += list(annotations.loc[annotations['fasta'] == fasta].index)

    # filter scaffolds
    if scaffolds is not None:
        for scaffold in scaffolds:
            genes_to_keep += list(annotations.loc[annotations['scaffold'] == scaffold].index)

    # filter genes
    if genes is not None:
        genes_to_keep += genes

    if (identifiers is not None) or (categories is not None):
        gene_to_ids = dict()
        for i, row in annotations.iterrows():
            gene_to_ids[i] = set(get_ids_from_row(row))
        # get genes with ids
        if identifiers is not None:
            identifiers = set(identifiers)
            for gene, ids in gene_to_ids.items():
                if len(set(ids) & set(identifiers)) > 0:
                    genes_to_keep.append(gene)

        # get genes from distillate categories
        if categories is not None:
            db_locs = get_database_locs()
            genome_summary_form = pd.read_csv(db_locs['genome_summary_form'], sep='\t')
            for level in ['module', 'sheet', 'header', 'subheader']:
                for category, frame in genome_summary_form.loc[~pd.isna(genome_summary_form['level'])].groupby(level):
                    if category in categories:
                        for gene, ids in gene_to_ids.items():
                            if len(ids & set(frame['gene_id'])) > 0:
                                genes_to_keep.append(gene)
                                break

    # make output
    output_fasta_generator = (i for i in read_sequence(input_fasta, format='fasta')
                              if i.metadata['id'] in set(genes_to_keep))
    write_sequence(output_fasta_generator, format='fasta', into=output_fasta)
