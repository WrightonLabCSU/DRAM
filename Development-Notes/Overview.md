# DRAM2 Nextflow Overview

Reed Woyda
Created June 24, 2024

This document serves as an overview of DRAM2. This will include motivations for DRAM2 (Rory) as well as motivations for DRAM2 (Nextflow). Additionally, this document will recap the status of both versions of DRAM2 and the progress which was made since development started for Reed Woyda. This document will cover most aspects of DRAM2 but will not contain how to run DRAM2 as this is currently within the GitHub repository main README file.

-------

## Table of Contents
- [DRAM2-Computational-Resource-Management](./DRAM2-Computational-Resource-Management.md)
- [DRAM2-Trees-Notes](./DRAM2-Trees-Notes.md)
- [DRAM2-Containers-Conda-Envs](./DRAM2-Containers-Conda-Envs.md)
- [DRAM2-Databases-Notes](./DRAM2-Databases-Notes.md)
- [DRAM2-Nextflow-Profiles](./DRAM2-Nextflow-Profiles.md)
- [DRAM2-RW-Daily-Logs](./DRAM2-RW-Daily-Logs.md)
- [How-To-Run-W2](./How-To-Run-W2.md)
- [How-To-Setup-W2-Riviera](./How-To-Setup-W2-Riviera.md)

------

## DRAM2 (Rory) Recap

A version of DRAM2 was initially started by Rory Flynn and the state of this pipeline will be described in this section. This version of DRAM2 was implemented almost entirely in python and was built from the ground up. The documentation for this version of DRAM2 is within the DRAM2-Archived repository within the WrightonLabCSU GitHub page. The documentation consists of the DRAM2-Archived repository and the associated ReadTheDocs documentation within the GitHub repository. Upon testing of DRAM2 in July/August 2023, to understand the current state of its functionality, it was determined that the DRAM2 functions for Call, Annotate, Distill, and Product were functional. It was assumed the Prepare Databases was functional but testing was not performed to further interrogate the status of this functionality. Other DRAM2 processes such as Traits/Adjectives, Trees, Neighborhoods, Strainer and others, were partially implemented but the status of these is unfinished.

Running this version of DRAM2 (Rory) is possible and can be used to call genes, annotate those genes and to generate both distill and product output. However, the state of the remaining functionality described in the help menu is unclear. More in-depth testing is needed to understand the current state of this version of DRAM2. Additionally, there do exist some bugs surrounding conda environment dependency versions. Also, the generation of a product seems to not include CAZY-annotated genes. I was unaware of this until I shared some results from a run and they noted that the expected presence of CAZY annotations, within the product, was missing. 

------

### DRAM1 Recap, DRAM2 (Rory) Recap and Nextflow Motivations

DRAM1 is currently in a failed pytest state. The choice was made in Fall 2023 to not support DRAM1 usage. While some users are still able to install and run DRAM1, novice users may encounter difficulties. Attempts to update and release version 1.5.0 after this decision were performed, to support the release, and inclusion into DRAM1, of the CAMPER database. GitHub issues were responded to indicating that DRAM1 is no longer supported and users should await the release of DRAM2 for annotation and metabolic characterization of FASTA files. 

Initially, the development of DRAM2 was going to proceed with the DRAM2 (Rory) version. However, due to reasons listed here, a development decision was made to implement DRAM2 in Nextflow. Reasons for DRAM2 implementation into Nextflow are: 1) Nextflow offers built-in support/functionality for the building of version-controlled and containerized bioinformatics pipelines, 2) DRAM, in general, has many processes which can be parallelized, something which Nextflow handles with ease, 3) The modulatiry of coding in Nextflow aligns with the goal of making separate workflows for the various processes DRAM2 performs. Additionally, as Nextflow is language-independent, recycling of existing DRAM2 (Rory) code is possible. 

Major reasons for re-working DRAM2 in Nextflow was to prevent/preempt two major issues with DRAM1: 1) dependency issues within Conda environments which require a technician to respond to failed tests and issues which require updating Conda environment conflicts and issues. 2) Preparing databases has resulted in a majority of DRAM1-related issues. This issue is mainly due to updates to the databases which are downloaded by the prepare_databases() function of DRAM1. Over time, databases change in size, format and version. Thus, a technician must respond to issues and updates to the included databases to ensure users may prepare these without error. Taken together, a solution to this problem was formulated: 1) Containerize dependencies within Singularity Containers, as well as providing Conda environments as a backup to users who are on HPCs without Singularity installed. These containers will be hosted either on the SyLabs library or via GLOBUS. 2) Provide pre-built databases to users. This results in the need for a dedicated support technician to download and format both the input databases, HMM or MMSeqs, and the descriptions database when updates are made to the included databases. While this is ideal for a user, having the ability to download from GLOBUS the pre-formatted databases and potentially an updated DRAM2, it does put a support burden on the hosting lab. Currently, there is no working Nextflow-based code for the downloading or preparing of the databases. DRAM2 (Nextflow) relies on previously built versions of the databases, and the descriptions database, which were built using DRAM1. A potential path forward would be to rely on DRAM2 (Rory) to build and format these databases and the descriptions database. However, there is some intial code to download the databases, and generate HMM and/or MMSeqs2 databases, but this is in an incomplete state.

------

## DRAM2 (Nextflow) Recap

Upon the decision to implement DRAM2 in Nextflow, a GitHub repository was created and an initial goal was set to rebuild the DRAM2 base functionality: Call, Annotate and Distill. This initial goal would be a proof-of-concept for implementation in Nextflow and because of this, a decision was made to rely on the DRAM2 (Rory) formatted databases and SQL descriptions database. This will require either: 1) Using the existing DRAM2 (Rory) code to prepare databases in the future, or 2) reusing the previous DRAM2 code for preparing databases or rewriting this code in a Nextflow-based process. 

------

### GitHub and ReadTheDocs

The DRAM2 GitHub has 3 types branches: 1) main - this branch is to push final working updates to. 2) dev - this is the branch which is actively developed on and is always the most up-to-date. 3) viz branches - these are working branches for the visualization aspects of dram. This documentation will mention work done in the visualization branches but more comprehensive documentation can be found within those branches.

------

### Documentation for DRAM2

The documentation for DRAM2 (Nextflow) is within the GitHub repository directory `Development-Notes/`. Within this directory there is this document, `Overview.md`, as well as a documents for `Daily-Development-Notes` - a long day-by-day account of DRAM2 development, `DRAM2-Trees` - a description of DRAM2 Trees process, the status of this process and the path forward, `DRAM2-Databases` - a description of the big picture of how DRAM2 will provide pre-built databases, why we want to provide these databases in this format, and an overview of the various plans forward. 

------

### Locations on W2 Server

Pulled GitHub repository location and development location:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow`

This directory contains various test data within `test_data/`.

Pulled Repository:
`DRAM2-NF`

Database location:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/databases`

Backup database location:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-database-backup-06252024`

The current Singularity container, `DRAM2-Nextflow-Main-Container-March252024-V4.sif`, is within the containers directory. 

Within this directory there are also various results directories and the databases directory. The databases directory is within the .gitignore. This contains the working versions of the DRAM2 databases - individual annotation databases and the SQL descriptions database.

Central location for Singularity images - there is a location where some containers are housed and the long-term plan would be to set the location within the config files, for the containers, to be here:
`/home/opt/singularity_containers/`
Note: This is where COMET Singularity images reside.

Currently there is a backup of all Singularity images for DRAM2 located at:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-singularity-container-backups`

------

### DRAM2 (Nextflow) Process overview

#### DRAM2 Modules, Assets, Containers and Environments

Modules are individual processes which are invoked through the pipeline. For example, `CALL_GENES()` is a process which takes in a FASTA file and uses Prodigal to identify coding and non-coding sequences. Modules are organized, by broad DRAM2 functionality, within the `modules` directory within the GitHub repository. Each individual module directory may contain numerous `.nf` files. These are the Nextflow process definitions for each individual process within the DRAM2 pipeline. These are set within lines 34 to 110 in DRAM2.nf.

Assets are various scripts and configuration files needed to run DRAM2. This contains both Conda and Singularity Profile definitions - both for SLURM-based DRAM2 runs and non-SLURM-based DRAM2 runs. As of now, this is a holding place for all of the individual python, R, Perl, etc., scripts which some processes rely on. These individual files are not ideal the the future goal would be to organize these together into a toolkit which will be placed here and symlinked to each individual process (as it is done now but this is messy and needs to be cleaned up). These are given variable names within the nextflow.config file and are put into channels, in the setup code of DRAM2.nf, based on the users command-line options. Assets also includes the internal directory where DRAM2 Databases is build built. This will be covered more in-depth in the documentation for DRAM2 Databases. Additionally, there are directories for both viz (visualization) and trees (needed input files). 

Containers directory holds the various Singularity containers which DRAM2 relies on. Currently, DRAM2 utilizes a single container, `DRAM2-Nextflow-Main-Container-March252024-V4.sif`. This is located on W2 at: `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/containers` (along with older versions). 

------

#### Main DRAM2.nf script and future plans

##### Pre-checks and file loading

Within the main Nextflow script, DRAM2.nf, there are an initial ~1,000 lines which check for various inputs and command line arguments. 

Lines 1 to 20 are a manifest which hold information about the DRAM2 pipeline and has developer information and documentation locations. Line 30 enables DSL2. 

Lines 34-100 load the various modules and the paths to the modules are assumed to stem from the pulled DRAM2 repository and may need to be changed if a user is going to install this system-wide. For example, DRAM2 on the W2 server would ideally hold the modules and Singularity containers in a central location in /opt/. 

Lines 114 - 877 are used to establish proper inputs from the users and to ensure the arguments the user provides are accurate and expected.

Lines 886 - 949 are used to set up various command-line outputs which are used to inform the user of the arguments and potions they used. There are multiple of these and based on the users command-line arguments and options, the correct one will be displayed. For example, there are different menus for call, annotate and distill. 

Lines 952 - 1440 are the main workflow for the DRAM2 Nextflow pipeline. This information is detailed within the sections here.

Lines 1446 - 1493 hold the various versions, used within the Singularity containers, and can be invoked by the command line argument `--v` or `--version`. These lines are followed by the various help menus. 

The main workflow, at lines 952 - 1440, contains all current code for running call, annotate, distill and product. While these are all included in a single workflow, the eventual goal is to separate the various DRAM2 processes (Call, Annotate, Distill, Product, Trees, etc.) into separate workflows. This will highly modularize DRAM2. However, now each DRAM2 process is invoked via the command-line options for call, annotate and distill. 

------

##### Renaming and calling genes

DRAM2 operates under the assumption that the user is providing individual FASTA files which are denoted as `samples`. It is also assumed that the headers of each individual FASTA file are unique for that given file. If this is not the case, the user may use the `--rename` option to rename the sample FASTA file headers with a unique prefix based on the basename of the sample provided. This name is obtained from the input FASTA file filename. This is performed using BBTools rename.sh script.

Calling genes is essential as coding sequences must be identified before they can be annotated. Prodigal is used for this step. The setup-code ensures the user has provided a path to input FASTA files. Prodigal only runs on 1 CPU so this is a step which can be highly parallelized. After gene calling, all of the individual output channels, for each individual sample, are combined and sent to QUAST, to get genome stats, and to rRNA scan and a tRNA scan processes to obtain present rRNA and tRNA genes. Both tRNA and rRNA outputs are collected and formatted in the respective `_COLLECT` processes. 

------

##### Annotation

Annotation starts with an initial check to see if the user had previous annotations they want to merge. This is a standalone process which takes in a directory with a variable number of DRAM2 (Nextflow) annotation (TSV) files. This is distinct from add_annotations().

For general annotation, an empty channel is created, `formattedOutputChannels`, which will hold a set of formatted output annotations for each sample, for each database. 

Gene locations are presented in the output annotations (TSV) file. To generate these locations, for each sample and for each gene identified, we invoke the GENE_LOCS() process. This process will output a formatted TSV file containing each called gene and its location within a given contig. 

Currently, there are 13 databases the user may choose to annotate with. Using the appropriate command-line options, the user may specify these in various manners. Based on the users choice, flags for each database are set and subsequently used to invoke the various databases to be annotated against. Each database may have an HMM formatted database (HMMR/hmmscan) or an MMSeqs2 formatted database, or both. If statements are used to check the database annotation flags and then proceed to query the database and format the output. For MMSeqs2 database searches, the input FASTA file is used to make and index, MMSEQS_INDEX(), which is then placed into a channel, combined with the FASTA file, which is then used by the process to query the databases. In the case of KEGG, which has an MMSeqs2 database, the index channel and FASTA channel are combined and are used to invoke MMSEQS_SEARCH_KEGG() which runs hmmscan. Next, SQL_KEGG() is invoked to generate a formatted output TSV file which will have the KEGG annotations along with the database descriptions, associated with those annotations, using the SQL descriptions database. This process is similar for HMM database searches but there are additional parsing and formatting steps. Once all annotations have been performed, the various formatted output channels are combined into one, `collected_formatted_hits`, and sent to the process COMBINE_ANNOTATIONS() to be combined into a single formatted TSV file. Additionally, depending on the users command-line options and arguments, the processes for adding bin (MAG) quality and bin (MAG) taxonomy information are performed. The user may also use a parameter to add annotations from a previous run into the current one. A check is done for this and if invoked, the ADD_ANNOTATIONS() process occurs. Lastly, a process called COUNT_ANNOTATIONS() is invoked to count unique annotations across the formatted annotations (TSV) file. 

Grouped in with annotation is the ability to output GFF and GENBANK annotation files. The user may specify these on the command-line when running DRAM2 (Nextflow). This will invoke the GENERATE_GFF_GENBANK() process.

Note on reverse best hits (RBH): RBH is only implemented for KEGG however, this portion of the code, within the KEGG annotation process, is commented out as it adds a significant amount of time to run KEGG annotations. This code can be uncommented and used as-is. It is suggested to optimize this process before a release.

------

##### Distill

Distill in DRAM2 is the act of reducing an annotations (TSV) file down to annotations of interest. The Wrighton Lab has, in collaboration with other researchers and research groups, to provide expertly-curated distill sheets. These sheets contain gene annotations and other notes and descriptions provided by the expert curators. Distill, like annotation databases, has many options for which distill sheets to use. The user may provide these various sheets/collections of sheets on the command line. A process is initially invoked, `COMBINE_DISTILL()` which simply collects the sheets the user requested and then subsequently the main `DISTILL()` process is invoked to generate a distill output. 

The distill output is a multi-sheet XLSX document. The idea is that each desired distill sheet requested will result in a single output sheet in the XLSX document. In addition to these individual output sheets, there is a Genome_Stats sheet, the first output sheet, which lists various genome stats pertaining to each sample provided. This may include rRNA and tRNA information (counts) if the user either provides these or if these were ran previously in the pipeline. Additionally, this may include bin (MAG) taxonomy and bin (MAG) quality. Lastly, if the user did run rRNA and tRNA scans, or provided the file separately with the required command-line options, there will be individual output sheets for each of these in the XLSX document. 

While distill functions correctly it is inefficient and needs optimization. As of now, it relies on brute-force parsing and matching of the input distill sheets and the annotations TSV file.

------

##### Product

There exists placeholder code for running a product at line 1419. This, at times, has been functional, based on the work done in the viz branches, but has been uncommented until tested and verified.

------

##### Trees

A more in-depth description of DRAM2 Trees and how to build the trees, is in the documentation file, DRAM2-Trees-Notes.md.

DRAM2 Trees functionality is incorporated into the main DRAM2.nf. Trees is designed to be run be default, unless specified by the user. DRAM2 Trees aims to reconcile annotations which are not able to be reconciled with the databases provided. 

Integration with Distill: While DRAM2 trees works to insert a new column into the annotations TSV file, work needs to be done for distill to recognize these updated annotations and proceed accordingly.

Note on running DRAM2: It is advised, until DRAM2 Trees is complete, to use the command-line option `--no_trees` to skip DRAM2 Trees.

