<p align="center">
  <img src="assets/images/DRAM2_large.png" width="600" height="600">
</p>

<h1 align="center">
  <strong> ⚠️ DRAM2 is currently under active development and usage is at your own risk. ⚠️ </strong></br>
</h1>


DRAM2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic assembled genomes. DRAM2 annotates MAGs using [KEGG](https://www.kegg.jp/) (if provided by the user), [UniRef90](https://www.uniprot.org/), [PFAM](https://pfam.xfam.org/), [dbCAN](http://bcb.unl.edu/dbCAN2/), [RefSeq viral](https://www.ncbi.nlm.nih.gov/genome/viruses/), [VOGDB](http://vogdb.org/) and the [MEROPS](https://www.ebi.ac.uk/merops/) peptidase database as well as custom user databases. DRAM is run in two stages. First an annotation step to assign database identifiers to gene, and then a distill step to curate these annotations into useful functional categories. DRAM2 was implemented in [Nextflow](https://www.nextflow.io/) due to its innate scalability on HPCs and containerization, ensuring rigorous reproducibility and version control, thus making it ideally suited for high-performance computing environments. 

For more detail on DRAM and how DRAM works please see our [paper](https://academic.oup.com/nar/article/48/16/8883/5884738) as well as the [wiki](https://github.com/WrightonLabCSU/DRAM/wiki).

### DRAM2 Development Note

The DRAM development team is actively working on DRAM2. We do not anticipate adding any additional functionality to DRAM, i.e. DRAM1.
- Future updates will include:
  - Both support for Nextflow + Conda and Nextflow + Singularity (Note: Singularity is not well-supported for MAC.)
  - Pre-formatted 

 
----------

<a name="install"></a>
### Installation
1) Clone the DRAM2 GitHub Repository
2) Download Singularity container and pre-formatted databases
2) [Install Nextflow >= v23.04.2.5870](https://www.nextflow.io/docs/latest/getstarted.html)
3) [Install Singularity >= v3.7.0](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (to pull Singualrity images from SyLabs).

#### General Instructions:
```
git clone https://github.com/WrightonLabCSU/DRAM2.git
cd COMET
./pull_singularity_containers.py
./pull_databases_full.py
```

#### Note for use on HPC:
- DRAM2 has SLURM auto submission integrated into the pipeline.
  - To use this feature, ensure you have [SLURM](https://slurm.schedmd.com/quickstart_admin.html) on your cluster.
- **If you do not have SLURM and do not want to use SLURM, use the provided alternative config file: `nextflow-No-SLURM.config`.**
  - **To use this config, you need to add the following to your command: `-c nextflow-No-SLURM.config`.**

---------

### DRAM2 Databases
DRAM2 databases, unlike DRAM1 databases, will be pre-formatted and hosted online. Users of DRAM2 will need to 1) decide which databases suits their needs and 2) download DRAM2 databases via the provided `pull_databases_*.py scripts. However, these databases can be quite large and it is therefore important to look through the options below.

These databases rely on an SQL database of database descriptions which is provided in 3 different sizes based on ther user's needs.

#### All databases DRAM2 accommodates:
- [KEGG](https://www.genome.jp/kegg/)
    - Kyoto Encyclopedia of Genes and Genomes.
    - (140G)
- [dbCAN](http://bcb.unl.edu/dbCAN2/)
    - A database for automated carbohydrate-active enzyme annotation.
    - (202M)
- [Kofam](https://www.genome.jp/tools/kofamkoala/)
    - Customized HMM database of KEGG Orthologs (KOs).
    - (14G)
- [MEROPS](https://www.ebi.ac.uk/merops/)
    - A database of proteolytic enzymes and their substrates.
    - (3.6G)
- [Viral](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239)
    - RefSeq viral database.
    - (1.6G)
- [CAMPER](https://github.com/WrightonLabCSU/CAMPER)
    - Curated Annotations for Microbial (Poly)phenol Enzymes and Reactions.
    - (846M)
- [CANT-HYD](https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation)
    - Curated database of phylogeny-derived hidden markov models for annotation of marker genes involved in hydrocarbon degradation.
    - (877M)
- [FeGenie](https://github.com/Arkadiy-Garber/FeGenie)
    - Placeholder description.
    - (6.6M)
- [Sulfur](url_to_sulfur_placeholder)
    - Custom Sulfur database.
    - (1.7M)
- [Methyl](url_to_methyl_placeholder)
    - Hidden Markov models (HMMs) based on genes related to iron acquisition, storage, and reduction/oxidation in Bacteria and Archaea.
    - (52K)
- [UniRef](https://www.uniprot.org/help/uniref)
    - A comprehensive and non-redundant database of protein sequences.
    - (477G)
- [Pfam](https://pfam.xfam.org/)
    - A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).
    - (8.8G)
- [VOGDB](url_to_vogdb_placeholder)
    - Placeholder description.
    - (4.5G)

#### DRAM2 Database (sets)
**Big set**
- Includes all medium databases + UniRef 
- Excludes KEGG*
- 
`./pull_databases_full.py`
OR
Follow these instructions to pull manually via [GLOBUS](https://www.globus.org/).

**Routine Set**
- Includes: dbCAN, Kofam, MEROPS, Viral, CAMPER, CANT-HYD, FeGenie, Sulfur, Methyl, Pfam, VOGDB
- Excludes KEGG*
- Excludes UniRef
  
`./pull_databases_routine.py`
OR
Follow these instructions to pull manually via [GLOBUS](https://www.globus.org/).

**Minimal Set** 
- Includes: <TBD>

`./pull_databases_minimal.py`
OR
Follow these instructions to pull manually via [GLOBUS](https://www.globus.org/).

##### KEGG download and format
A [subscription](https://www.kegg.jp/kegg/download/) is required to download the [kegg](https://www.genome.jp/kegg/) database.
If you have a subscription and would like to use KEGG you will need to 1) download the KEGG <blank> file, 2) run the `format-KEGG-DRAM2.py script` to format the KEGG database for DRAM2 and 3) place the formatted database files in the appropriate location.

[Download](https://www.kegg.jp/kegg/download/) the amino acid FASTA file (`*.pep`) downloaded from KEGG. This can be any of the gene fasta files that are provided by the KEGG FTP server or a concatenated version of them. 

Run the `format-KEGG-DRAM2.sh` script:
`./format-KEGG-DRAM2.py <path-to-downloaded-kegg.pep>`

--------

#### DRAM2 Descriptions Database
**Big set**
- Includes Uniref and KEGG
`./pull_descriptions_full.py`
OR
Follow these instructions to pull manually via (GLOBUS)[https://www.globus.org/].

**Routine Set**
- Excludes Uniref and KEGG
`./pull_descriptions_routine.py`
OR
Follow these instructions to pull manually via (GLOBUS)[https://www.globus.org/].

--------

## System Requirements
**EDIT**


*NOTE:* Setting up DRAM can take a long time (up to 5 hours) and uses a large amount of memory (512 gb) by default. To
use less memory you can use the `--skip_uniref` flag which will reduce memory usage to ~64 gb if you do not provide KEGG
 Genes and 128 gb if you do. Depending on the number of processors which you tell  it to use (using the `--threads`
argument) and the speed of your internet connection. On a less than 5 year old server with 10 processors it takes about
 2 hours to process the data when databases do not need to be downloaded.



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


