import pandas as pd
from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.summarize_vgfs import get_ids_from_row

# TODO: filter by taxonomic level, completeness, contamination
# TODO: filter scaffolds file, gff file
# TODO: add negate, aka pull not from given list
# TODO: update get_ids_from_row to get PFAM ids
# TODO: pull sequences based on distillate categories

def pull_sequences(input_annotations, input_fasta, output_fasta, fastas=None, scaffolds=None, genes=None,
                   identifiers=None):
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

    # get genes with ids
    if identifiers is not None:
        identifiers = set(identifiers)
        for i, row in annotations.iterrows():
            row_ids = get_ids_from_row(row)
            if len(set(row_ids) & set(identifiers)) > 0:
                genes_to_keep.append(i)

    # make output
    output_fasta_generator = (i for i in read_sequence(input_fasta, format='fasta') if i.metadata['id'] in genes_to_keep)
    write_sequence(output_fasta_generator, format='fasta', into=output_fasta)
