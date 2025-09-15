#!/usr/bin/env python

import os
import pandas as pd
from pathlib import Path

FASTA_COLUMN = os.getenv('FASTA_COLUMN')

df = pd.read_csv(Path("combined_trna_scan.tsv").resolve(), sep='\t')

gene_col="gene_id"
df = df[[FASTA_COLUMN, gene_col]].copy()
counts = pd.crosstab(df[gene_col], df[FASTA_COLUMN], dropna=False).astype("int64")
counts.columns.name = None
counts = counts.reset_index()
counts["gene_description"] = counts[gene_col] + " tRNA"
counts["category"] = "tRNA"
counts["topic_ecosystem"] = "tRNA"
counts["subcategory"] = ""

counts.to_csv("collected_trnas.tsv", sep="\t", index=False)
