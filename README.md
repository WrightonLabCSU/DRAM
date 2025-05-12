# DRAM v2

<p align="center">
  <img src="assets/images/DRAM2_large.png" width="600" height="600" alt="DRAM v2 logo">
</p>

## ⚠️ DRAM v2 is currently under active development and usage is at your own risk. ⚠️

DRAM v2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic and genomic assembled data (e.g. scaffolds or contigs) or called genes (e.g. nuclotide or amino acid format). DRAM annotates MAGs using [KEGG](https://www.kegg.jp/) (if provided by the user), [UniRef90](https://www.uniprot.org/), [PFAM](https://pfam.xfam.org/), [dbCAN](http://bcb.unl.edu/dbCAN2/), [RefSeq viral](https://www.ncbi.nlm.nih.gov/genome/viruses/), [VOGDB](http://vogdb.org/) and the [MEROPS](https://www.ebi.ac.uk/merops/) peptidase database as well as custom user databases.

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

- [Installation Guide](docs/installation.md)
- [Usage Examples](#example-usage)
- [Command Line Options](#command-line-options)
- [Nextflow Tips and Tricks](#nextflow-tips-and-tricks)
- [Resource Management](#resource-management)
- [Citation](#citing-dram)

## Example Usage

DRAM apps Call, Annotate and Distill can all be run at once or alternatively, each app can be run individually. Here are some common usage examples:

1) **Rename fasta headers based on input sample file names:**
```bash
nextflow run DRAM --rename --input_fasta <path/to/fasta/directory/>
```

2) **Call genes using input fastas (use --rename to rename FASTA headers):**
```bash
nextflow run DRAM --call --rename --input_fasta <path/to/fasta/directory/>
```

3) **Annotate called genes using input called genes and the KOFAM database:**
```bash
nextflow run DRAM --annotate --input_genes <path/to/called/genes/directory> --use_kofam
```

4) **Annotate called genes using input fasta files and the KOFAM database:**
```bash
nextflow run DRAM --annotate --input_fasta <path/to/called/genes/directory> --use_kofam
```

5) **Merge various existing annotations files together (Must be generated using DRAM):**
```bash
nextflow run DRAM --merge_annotations <path/to/directory/with/multiple/annotation/TSV/files>
```

6) **Distill using input annotations:**
```bash
nextflow run DRAM --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>
```

7) **Complete workflow example:**
```bash
nextflow run -bg DRAM \
  --input_fasta ../test_data/DRAM_test_data/ \
  --outdir DRAM-test-data-call-annotate-distill \
  --threads 8 \
  --call --rename --annotate \
  --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper \
  --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur \
  --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom test-data/custom-test-distilalte.tsv \
  -profile conda_slurm --slurm_node main -with-report -with-trace -with-timeline
```

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