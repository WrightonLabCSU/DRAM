# DRAM2 Nextflow Overview

Reed Woyda
Created June 24, 2024

This document serves as an overview of DRAM2. This will include motivations for DRAM2 (Rory) as well as motivations for DRAM2 (Nextflow). Additionally, this document will recap the status of both versions of DRAM2 and the progress which was made since development started for Reed Woyda. This document will cover most aspects of DRAM2 but will not contain how to run DRAM2 as this is currently within the GitHub repository main README file.


## DRAM2 (Rory) Recap

A version of DRAM2 was initially started by Rory Flynn and the state of this pipeline will be described here. This version of DRAM2 was implemented almost entirely in python and was built from the ground up. The documentation for this version of DRAM2 is within the DRAM2-Archived repository within the WrightonLabCSU GitHub page. The documentation consists of the DRAM2-Archived repository and the associated ReadTheDocs documentation within the GitHub repository. Upon testing of DRAM2 in July/August 2023, to understand the current state of its functionality, it was determined that the DRAM2 functions for Call, Annotate, Distill, and Product were functional. It was assumed the Prepare Databases was functional but testing was not performed to further interrogate the status of this functionality. Other DRAM2 processes such as Traits/Adjectives, Trees, Neighborhoods, Strainer and others, were partially implemented but the status of these is unfinished. 

### DRAM2 (Rory) Recap and Nextflow Motivations

Initially, the development of DRAM2 was going to proceed with the DRAM2 (Rory) version. However, due to reasons listed here, a development decision was made to implement DRAM2 in Nextflow. Reasons for DRAM2 implementation into Nextflow are: 1) Nextflow offers built-in support/functionality for the building of version-controlled and containerized bioinformatics pipelines, 2) DRAM, in general, has many processes which can be parallelized, something which Nextflow handles with ease, 3) The modulatiry of coding in Nextflow aligns with the goal of making separate workflows for the various processes DRAM2 performs. Additionally, as Nextflow is language-independent, recycling of existing DRAM2 (Rory) code is possible. 

## DRAM2 (Nextflow) Recap

Upon the decision to implement DRAM2 in Nextflow, a GitHub repository was created and the initial goal was to rebuild DRAM2 base functionality: Call, Annotate and Distill. This initial goal would be a proof-of-concept for implementation in Nextflow and because of this, a decision was made to rely on the DRAM2 (Rory) formatted databases and SQL descriptions database. This will require either: 1) Using the existing DRAM2 (Rory) code to prepare databases in the future, or 2) reusing the previous DRAM2 code for preparing databases or rewriting this code in a Nextflow-based process. 

### GitHub and ReadTheDocs

The DRAM2 GitHub has 3 main branches: 1) main - this branch is to push final working updates to. 2) dev - this is the branch which is actively developed on and is always the most up-to-date. 3) viz branches - these are working branches for the visualization aspects of dram. This documentation will mention work done in the visualization branches but more comprehensive documentation can be found within those branches.

### DRAM2 (Nextflow) Process overview

#### DRAM2 Modules, Assets, Containers and Environments

Modules are individual processes which are invoked through the pipeline. For example, `CALL_GENES()` is a process which takes in a FASTA file and uses Prodigal to identify coding and non-coding sequences. Modules are organized, by broad DRAM2 functionality, within the `modules` directory within the GitHub repository. Each individual module directory may contain numerous `.nf` files. These are the Nextflow process definitions for each individual process within the DRAM2 pipeline.

Assets are various scripts and configuration files needed to run DRAM2. This contains both Conda and Singularity Profile definitions - both for SLURM-based DRAM2 runs and non-SLURM-based DRAM2 runs. As of now, this is a holding place for all of the individual python, R, Perl, etc., scripts which some processes rely on. These individual files are not ideal the the future goal would be to organize these together into a toolkit which will be placed here and symlinked to each individual process (as it is done now but this is messy and needs to be cleaned up). Assets also includes the internal directory where DRAM2 Databases is build built. This will be covered more in-depth in the documentation for DRAM2 Databases. Additionally there are directories for both viz (visualization) and trees (needed input files).

Containers directory holds the various Singularity containers which DRAM2 relies on 

#### Main Workflow and future plans

The main workflow, at lines 952 - 1440, contains all current code for running call, annotate, distill and product. While these are all included in a single workflow, the eventual goal is to separate the various DRAM2 processes (Call, Annotate, Distill, Product, Trees, etc.) into separate workflows. This will highly modularize DRAM2.

##### Renaming and calling genes

DRAM2 operates under the assumption that the user is providing individual FASTA files which are denoted as `samples`. It is also assumed that the headers of each individual FASTA file are unique for that given file. If this is not the case, the user may use the `--rename` option to rename the sample FASTA file headers with a unique prefix based on the basename of the sample provided. This name is obtained from the input FASTA file filename. This is performed using BBTools rename.sh script.

Calling genes is essential as coding sequences must be identified before they can be annotated. Prodigal is used for this step


#### Pre-checks and file loading

Within the main Nextflow script, DRAM2.nf, there are an initial ~1,000 lines which check for various inputs and command line arguments. 

Lines 1 to 20 are a manifest which hold information about the DRAM2 pipeline and has developer information and documentation locations. Line 30 enables DSL2. 

Lines 34-100 load the various modules and the paths to the modules are assumed to stem from the pulled DRAM2 repository and may need to be changed if a user is going to install this system-wide. For example, DRAM2 on the W2 server would ideally hold the modules and Singularity containers in a central location in /opt/. 


Lines 114 - 877 are used to establish proper inputs from the users and to ensure the arguments the user provides are accurate and expected.

Lines 886 - 949 are used to set up various command-line outputs which are used to inform the user of the arguments and potions they used. There are multiple of these and based on the users command-line arguments and options, the correct one will be displayed. For example, there are different menus for call, annotate and distill. 

Lines 952 - 1440 are the main workflow for the DRAM2 Nextflow pipeline.

Lines 1446 - 1493 hold the various versions, used within the Singularity containers, and can be invoked by the command line argument `--v` or `--version`. These lines are followed by the various help menus. 






### Documentation for DRAM2

The documentation for DRAM2 (Nextflow) is within the GitHub repository directory 'Development-Notes`. Within this directory there is this document, `Overview.md`, as well as a documents for `Daily-Development-Notes` - a long day-by-day account of DRAM2 development, `DRAM2-Trees` - a description of DRAM2 Trees process, the status of this process and the path forward, `DRAM2-Databases` - a description of the big picture of how DRAM2 will provide pre-built databases, why we want to provide these databases in this format, and an overview of the various plans forward. 