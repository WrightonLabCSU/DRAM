#!/usr/bin/env python

import os
import pandas as pd
from pathlib import Path

FASTA_COLUMN = os.getenv('FASTA_COLUMN', 'input_fasta')

df = pd.read_csv(Path("raw_trna_scan.tsv").resolve(), sep='\t')

gene_col="gene_id"
df = df[[FASTA_COLUMN, gene_col]].copy()
counts = pd.crosstab(df[gene_col], df[FASTA_COLUMN], dropna=False).astype("int64")
counts.columns.name = None
counts = counts.reset_index()
counts["AA_type"] = counts[gene_col].str.split(" ").str[0]
counts["gene_description"] = counts[gene_col] + " tRNA"
counts["category"] = "tRNA"
counts["topic_ecosystem"] = "tRNA"
counts["subcategory"] = ""

counts.to_csv("collected_trnas.tsv", sep="\t", index=False)
