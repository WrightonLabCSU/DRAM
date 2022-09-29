# DRAM
[![CircleCI](https://circleci.com/gh/WrightonLabCSU/DRAM/tree/master.svg?style=svg)](https://circleci.com/gh/WrightonLabCSU/DRAM/tree/master)

DRAM (Distilled and Refined Annotation of Metabolism) is a tool for annotating metagenomic assembled genomes and [VirSorter](https://github.com/simroux/VirSorter) identified viral contigs. DRAM annotates MAGs and viral contigs using [KEGG](https://www.kegg.jp/) (if provided by the user), [UniRef90](https://www.uniprot.org/), [PFAM](https://pfam.xfam.org/), [dbCAN](http://bcb.unl.edu/dbCAN2/), [RefSeq viral](https://www.ncbi.nlm.nih.gov/genome/viruses/), [VOGDB](http://vogdb.org/) and the [MEROPS](https://www.ebi.ac.uk/merops/) peptidase database as well as custom user databases. DRAM is run in two stages. First an annotation step to assign database identifiers to gene, and then a distill step to curate these annotations into useful functional categories. Additionally, viral contigs are further analyzed during to identify potential AMGs. This is done via assigning an auxiliary score and flags representing the confidence that a gene is both metabolic and viral.

For more detail on DRAM and how DRAM works please see our [paper](https://academic.oup.com/nar/article/48/16/8883/5884738) as well as the [wiki](https://github.com/shafferm/DRAM/wiki).
For information on how DRAM is changing, please read the [release note](https://github.com/WrightonLabCSU/DRAM/releases/latest)

## Getting Started Part 1: Installation

**NOTE** If you already have an old release of DRAM installed and just want to upgrade, then please read the set-up step before you remove your old environment.

To install DRAM you also must install some dependencies. The easiest way to install both DRAM and its dependencies is to use [conda](https://docs.conda.io/en/latest/miniconda.html),  but you can also use manual instructions, or if you are an adventurer you can install a release candidate from this repository .
   
### Conda Installation

Install DRAM into a new [conda](https://docs.conda.io/en/latest/) environment using the provided
environment.yaml file.
```bash
wget https://raw.githubusercontent.com/shafferm/DRAM/master/environment.yaml
conda env create -f environment.yaml -n DRAM
```
If this installation method is used, then all further steps should be run inside the newly created DRAM environment, or with the full path to the executable, use `which` with the environment active to find these, the eg. `which DRAM.py`. This environment can be activated using this command:
```bash
conda activate DRAM
```

You have now installed DRAM, and are ready to set up the databases.

### Manual Installation

If you do not install via a conda environment, then the dependencies [pandas](https://pandas.pydata.org/), [networkx](https://networkx.github.io/), [scikit-bio](http://scikit-bio.org/), [prodigal](https://github.com/hyattpd/Prodigal), [mmseqs2](https://github.com/soedinglab/mmseqs2), [hmmer](http://hmmer.org/) and [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) need to be installed manually. Then you can install DRAM using pip:
```bash
pip install DRAM-bio
```

You have now installed DRAM, and are ready to set up the databases.

### Release Candidate Installation

The latest version of DRAM is often a release candidate, and these are not pushed to pypi, or Bioconda and so can't be installed with the methods above. You can tell if there is currently a release candidate by reading the [release notes](https://github.com/WrightonLabCSU/DRAM/releases).

To install a potentially unstable release candidate, follow the instructions below. Note the comments within the code sections as there is a context in which commands must be used.

```bash
# Clone the git repository and move into it
git clone https://github.com/WrightonLabCSU/DRAM.git
cd DRAM
# Install dependencies, this will also install a stable version of DRAM that will then be replaced.
conda env create --name my_dram_env -f environment.yaml
conda activate my_dram_env
# Install pip
conda install pip3
pip3 install ./
```

You have now installed DRAM, and are ready to set up the databases.


## Getting Started Part 2: Setup Databases

### I Want to Use an Already Setup Databases

If you already installed and set up a previous version of dram and want to use your old databases, then you can do it with two steps.

Activate your old DRAM environment, and save your old config:

```bash
conda activate my_old_env
DRAM-setup.py export_config > my_old_config.txt
```

Activate your new DRAM environment, and import your old databases

```bash
conda activate my_new_env
DRAM-setup.py import_config --config_loc  my_old_config.txt
```

### I have access to KEGG

Set up DRAM using the following command:

```bash
DRAM-setup.py prepare_databases --output_dir DRAM_data --kegg_loc kegg.pep
```

`kegg.pep` is the path to the amino acid FASTA file downloaded from KEGG. This can be any of the gene fasta files that are provided by the KEGG FTP server or a concatenated version of them. `DRAM_data` is the path  to the processed databases used by DRAM. If you already have any of the databases downloaded to your server and don't want to download them again then you can pass them to the `prepare_databases` command by use the `--{db_name}_loc` flags such as `--uniref_loc` and `--viral_loc`.

### I don't have access to KEGG

Not a problem. Then use this command:

```bash
DRAM-setup.py prepare_databases --output_dir DRAM_data
```

Similar to above you can still provide locations of databases you have already downloaded so you don't have to do it
again.

To test that your set up worked use the command `DRAM-setup.py print_config` and the location of all databases provided
will be shown as well as the presence of additional annotation information.

*NOTE:* Setting up DRAM can take a long time (up to 5 hours) and uses a large amount of memory (512 gb) by default. To
use less memory you can use the `--skip_uniref` flag which will reduce memory usage to ~64 gb if you do not provide KEGG
 Genes and 128 gb if you do. Depending on the number of processors which you tell  it to use (using the `--threads`
argument) and the speed of your internet connection. On a less than 5 year old server with 10 processors it takes about
 2 hours to process the data when databases do not need to be downloaded.

## Getting Started Part 3: Usage

Once DRAM is set up you are ready to annotate some MAGs. The following command will generate your full annotation:

```bash
DRAM.py annotate -i 'my_bins/*.fa' -o annotation
```

`my_bins` should be replaced with the path to a directory which contains all of your bins you would like to annotated and `.fa` should be replaced with the file extension used for your bins (i.e. `.fasta`, `.fna`, etc). If you only need to annotate a single genome (or an entire assembly) a direct path to a nucleotide fasta should be provided. Using 20 processors, DRAM.py takes about 17 hours to annotate ~80 MAGs of medium quality or higher from a mouse gut metagenome.

In the output `annotation` folder, there will be various files. `genes.faa` and `genes.fna` are fasta files with all genes called by prodigal with additional header information gained from the annotation as nucleotide and amino acid records respectively. `genes.gff` is a GFF3 with the same annotation information as well as gene locations. `scaffolds.fna` is a collection of all scaffolds/contigs given as input to `DRAM.py annotate` with added bin information in the headers. `annotations.tsv` is the most important output of the annotation. This includes all annotation information about every gene from all MAGs. Each line is a different gene and each column contains annotation information. `trnas.tsv` contains a summary of the tRNAs found in each MAG.

Then after your annotation is finished you can summarize these annotations with the following command:

```bash
DRAM.py distill -i annotation/annotations.tsv -o genome_summaries --trna_path annotation/trnas.tsv --rrna_path annotation/rrnas.tsv
```
This will generate the distillate and liquor files.

## System Requirements
DRAM has a large memory burden and is designed to be run on high performance computers. DRAM annotates against a large
variety of databases which must be processed and stored. Setting up DRAM with KEGG Genes and UniRef90 will take up ~500
GB of storage after processing and require ~512 GB of RAM while using KOfam and skipping UniRef90 will mean all
processed databases will take up ~30 GB of disk and will only use ~128 GB of RAM while processing. DRAM annotation
memory usage depends on the databases used. When annotating with UniRef90, around 220 GB of RAM is required. If the KEGG
gene database has been provided and UniRef90 is not used, then memory usage is around 100 GB of RAM. If KOfam is used to
annotate KEGG and UniRef90 is not used, then less than 50 GB of RAM is required. DRAM can be run with any number of
processors on a single node.

## Citing DRAM
The DRAM was published in Nucleic Acids Research in 2020 and is available [here](https://academic.oup.com/nar/article/48/16/8883/5884738). If
DRAM helps you out in your research, please cite it.


