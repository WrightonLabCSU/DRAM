# DRAM v2

<p align="center">
  <img src="assets/images/DRAM2_large.png" width="600" height="600" alt="DRAM v2 logo">
</p>

## ⚠️ DRAM v2 is currently under active development and usage is at your own risk. ⚠️

DRAM v2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic and genomic assembled data (e.g. scaffolds or contigs) or called genes (e.g. nuclotide or amino acid format). DRAM annotates MAGs using [KEGG](https://www.kegg.jp/) (if provided by the user), [UniRef90](https://www.uniprot.org/), [PFAM](https://pfam.xfam.org/), [dbCAN](http://bcb.unl.edu/dbCAN2/), [RefSeq viral](https://www.ncbi.nlm.nih.gov/genome/viruses/), [VOGDB](http://vogdb.org/) and the [MEROPS](https://www.ebi.ac.uk/merops/) peptidase database as well as custom user databases.

Viral catalogs from a typical geNomad → CheckV pipeline are also supported via **DRAM-v Phase 1 viral mode** (`--use_dramv true`): per-vMAG (per-scaffold) AMG flagging (M/K/E/A/P/T/F/B per the v1 conventions) and a viral-flavoured distillate, no VirSorter affi-contigs file required. See the "Viral mode" example below.

DRAM is run in four stages: 
1) Gene Calling Prodogal - genes are called on user provided scaffolds or contigs 
2) Gene Annotation - genes are annotated with a set of user defined databases 
3) Distillation - annotations are curated into functional categories
4) Product Generation - interactive visualizations of DRAM output are generated 

For more detail on DRAM and how DRAM v2 works please see our DRAM products:
- [DRAM version 1 publication](https://academic.oup.com/nar/article/48/16/8883/5884738)
- [DRAM in KBase publication](https://pubmed.ncbi.nlm.nih.gov/36857575/)
- [DRAM webinar](https://www.youtube.com/watch?v=-Ky2fz2vw2s)

## Quick Links

- [Docs](https://dramit.readthedocs.io/en/latest)
- [Installation Guide](https://dramit.readthedocs.io/en/latest/installation.html)
- [Usage Examples](https://dramit.readthedocs.io/en/latest/usage.html)
- [Parameter API]([#command-line-options](https://dramit.readthedocs.io/en/latest/params_doc.html))
- [Rules API]([#nextflow-tips-and-tricks](https://dramit.readthedocs.io/en/latest/rules_parser.html))
- [DRAM-v (Viral mode)](#dram-v-viral-mode) — per-sample / catalog launches, flag reference, `auxiliary_score`, geNomad adapter

## Example Usage

DRAM apps Call, Annotate and Distill can all be run at once or alternatively, each app can be run individually. Here are some common usage examples:

1) **Rename fasta headers based on input sample file names:**
```bash
nextflow run WrightonLabCSU/DRAM --rename --input_fasta <path/to/fasta/directory/>
```

2) **Call genes using input fastas (use --rename to rename FASTA headers):**
```bash
nextflow run WrightonLabCSU/DRAM --call --rename --input_fasta <path/to/fasta/directory/>
```

3) **Annotate called genes using input called genes and the KOFAM database:**
```bash
nextflow run WrightonLabCSU/DRAM --annotate --input_genes <path/to/called/genes/directory> --use_kofam
```

4) **Annotate called genes using input fasta files and the KOFAM database:**
```bash
nextflow run WrightonLabCSU/DRAM --annotate --input_fasta <path/to/called/genes/directory> --use_kofam
```

5) **Merge various existing annotations files together (Must be generated using DRAM):**
```bash
nextflow run WrightonLabCSU/DRAM --merge_annotations <path/to/directory/with/multiple/annotation/TSV/files>
```

6) **Distill using input annotations:**
```bash
nextflow run WrightonLabCSU/DRAM --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>
```

7) **Complete workflow example:**
```bash
nextflow run -bg WrightonLabCSU/DRAM \
  --input_fasta [DIRECTORY of fasta files] \
  --outdir [OUTPUT] \
  --rename --sum_ecos 'eng_sys,ag' \
  -profile singularity,full_mode
```

8) **Viral mode (DRAM-v) — AMG flags on geNomad+CheckV viral contigs.** See the dedicated [DRAM-v (Viral mode)](#dram-v-viral-mode) chapter below for the full guide (per-sample vs catalog, flag reference, `auxiliary_score`, and the geNomad adapter).

## DRAM-v (Viral mode)

DRAM-v adds three columns to the per-gene table and ships an AMG-filtered distillate. Designed for input from a geNomad → CheckV pipeline (no VirSorter affi-contigs file required). Two run modes — pick whichever matches your input:

- **Per-sample mode** — DRAM runs separately on each sample's filtered viral fasta (`<sample>_filtered.fna`). Used during development, smoke testing, or when you want per-sample annotations without going through clustering.
- **Catalog mode** — DRAM runs once on a clustered vOTU catalog whose contigs were renamed `<sample>_<orig>` upstream. Standard MIUViG / Sullivan-lab production path, recommended for any multi-sample analysis.

### What `--use_dramv` produces

Three new columns in `ANNOTATE/raw-annotations.tsv` and `ANNOTATE/annotations_with_flags.tsv`:

| Column | Type | Meaning |
| --- | --- | --- |
| `amg_flags` | str | concatenated single-letter flags in v1 order: `M K E V A P T F B`, plus a non-v1 `N` (see below) |
| `is_transposon` | bool | the gene's Pfam hits intersect a curated transposon set |
| `auxiliary_score` | int 1-5 | v1 flank-confidence score; lower = stronger viral context |

The `--amg_only` distillate (a single `AMG` sheet in `metabolism_summary.xlsx`, one count column per scaffold) keeps rows where `amg_flags` contains `M`, lacks `A`/`P`/`T`/`N`, **and** `auxiliary_score ≤ --max_auxiliary_score` (default 3). Set `--max_auxiliary_score 5` to disable the score filter.

Viral mode forces `groupby_column=scaffold` and skips QUAST + rRNA/tRNA collection — none of those make sense per-scaffold.

### Per-sample mode

```bash
nextflow run WrightonLabCSU/DRAM \
  --input_fasta path/to/filtered_fastas_dir \
  --fasta_fmt "*.fna" \
  --outdir results/dramv \
  --call --annotate --summarize \
  --use_kofam --use_dbcan --use_merops \
  --use_dramv true \
  --genomad_genes "path/to/genomad_genes/*_virus_genes.tsv" \
  -profile singularity
```

Each gene id matches between DRAM and geNomad (both call genes from scratch on the same per-sample fasta), so the geNomad → DRAM join hits via direct gene-id equality.

### Catalog mode

```bash
nextflow run WrightonLabCSU/DRAM \
  --input_fasta path/to/votu_catalog.fa \
  --outdir results/dramv \
  --call --annotate --summarize \
  --use_kofam --use_dbcan --use_merops \
  --use_dramv true \
  --genomad_genes "path/to/genomad_genes/*_virus_genes.tsv" \
  --genomad_filename_prefix true \
  -profile singularity
```

Catalog contigs are typically renamed `<sample>_<orig>` upstream so cross-sample names don't collide. Without `--genomad_filename_prefix true`, geNomad's per-sample gene ids (`<orig>_<num>`) won't match DRAM's catalog gene ids (`<sample>_<orig>_<num>`) and `auxiliary_score` collapses to 5 for every row — which the default `--max_auxiliary_score 3` filter would then drop entirely. With the flag on, each `<sample>_virus_genes.tsv` filename is parsed to derive the prefix and the join lines up.

### Flag reference

| Flag | Meaning |
| --- | --- |
| **`M`** | metabolism — gene matches a curated metabolic-gene set (KEGG/Pfam/CAZy/EC) |
| **`K`** | gene id appears in `bin/assets/amg_database.tsv`. `K` force-sets `M` per v1 |
| **`E`** | same as `K`, restricted to verified AMG rows |
| **`V`** | gene has a VOG hit whose VOGdb `FunctionalCategory` is `Xr` (replication) or `Xs` (structure). Requires `vog_annotations_latest.tsv`; `--use_dramv` auto-enables `--use_vog` |
| **`A`** | gene matches the curated cell-entry CAZy set |
| **`P`** | gene matches the curated viral peptidase MEROPS set |
| **`T`** | any gene on the same scaffold has `is_transposon=true` |
| **`F`** | gene is within `--amg_length_from_end` (default 5000) bp of either contig end |
| **`B`** | gene is part of a 3-consecutive-`M`-gene window on the scaffold |
| **`N`** *(not v1)* | gene id hits `amg_database.tsv` rows where `essential_viral_function=TRUE`, per [Martin et al. 2025](https://doi.org/10.1038/s41564-025-02095-4). Paper-cautioned genes (DsrC, QueC/QueF, folA/folB/folK, RNR, mazG, pur*, etc.) likely essential for viral processes rather than auxiliary metabolism. `N`-flagged rows are excluded from the strict-AMG distillate but stay in the full `raw-annotations.tsv` for review. `N` does **not** force `M`. |

### `auxiliary_score` reference

Verbatim port of v1's flank-confidence algorithm (`mag_annotator/annotate_vgfs.py:calculate_auxiliary_scores`, commit `6cd68f9`). Lower is better:

| Score | Condition |
| --- | --- |
| 1 | hallmark virus genes on **both** flanks |
| 2 | hallmark on one side, viral-like on the other |
| 3 | viral-like on both sides |
| 4 | hallmark/viral-like on at least one side, **or** self carries one |
| 5 | first/last on scaffold, **or** no viral context anywhere |

Then the `B`-flag downgrade fires: any gene with `B` in `amg_flags` and `auxiliary_score < 4` is bumped to 4 (a stretch of three metabolic genes is suspicious in a true viral region).

Without `--genomad_genes`, every gene scores 5 (v1 fallback for genomes without VirSorter input). Provide it to populate the score: geNomad's `virus_hallmark = 1` rows and marker suffix `.VV` / `.Vv` map to VirSorter hallmark (cat `0`); marker suffix `.vV` / `.vv` maps to viral-like (cat `1`).

The DRAM-v ↔ geNomad join is hybrid: direct gene-id equality first, position-overlap on the parent assembly contig as fallback. CheckV-trimmed provirus contigs (`<contig>|provirus_<start>_<end>`) are handled in both paths.

## Nextflow Tips and Tricks

The `-resume` option in Nextflow DSL2 allows you to efficiently manage and modify your workflow runs:

- **Adding databases to an existing run:**
  - Using `-resume` with your existing work directory lets you reuse called genes and existing annotations
  - Example: If you initially used `--use_kofam --use_dbcan`, you can add `--use_kegg --use_uniref` and only the new annotations will be computed

## Resource Management

DRAM leverages Nextflow's horizontal scaling capabilities to distribute computational tasks across multiple computing resources. You can customize resource allocation through the `nextflow.config` file:

- Modify "maxForks" parameters to control parallel execution
- Configure CPU and memory requirements per process
- Coming soon: "lite", "medium" and "heavy" modes for different computing environments

## Configuration

Every CLI option can be set in the `nextflow.config` file. For example:

```nextflow
params {
    use_uniref = true
    annotate = true
}
```

You can also use a custom config file:
```bash
nextflow run DRAM -c /path/to/custom_config.config
```

## Citing DRAM

If DRAM helps you in your research, please cite:
[DRAM publication in Nucleic Acids Research (2020)](https://academic.oup.com/nar/article/48/16/8883/5884738)
