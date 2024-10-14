<p align="center">
  <img src="assets/images/DRAM2_large.png" width="600" height="600" alt="DRAM v2 logo">
</p>

## ⚠️ DRAM v2 is currently under active development and usage is at your own risk. ⚠️


DRAM v2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic and genomic assembled data (e.g. scaffolds or contigs) or called genes (e.g. nuclotide or amino acid format). DRAM annotates MAGs using [KEGG](https://www.kegg.jp/) (if provided by the user), [UniRef90](https://www.uniprot.org/), [PFAM](https://pfam.xfam.org/), [dbCAN](http://bcb.unl.edu/dbCAN2/), [RefSeq viral](https://www.ncbi.nlm.nih.gov/genome/viruses/), [VOGDB](http://vogdb.org/) and the [MEROPS](https://www.ebi.ac.uk/merops/) peptidase database as well as custom user databases. DRAM is run in two stages. First an annotation step to assign database identifiers to gene, and then a distill step to curate these annotations into useful functional categories. DRAM was implemented in [Nextflow](https://www.nextflow.io/) due to its innate scalability on HPCs and containerization, ensuring rigorous reproducibility and version control, thus making it ideally suited for high-performance computing environments. 

DRAM is run in four stages: 
1) Gene calling - genes are called on user provided scaffolds or contigs 
2) Gene annotation - genes are annotated with a set of user defined databases 
3) Distillation - annotations are curated into functional categories
4) Product generation - interactive visualizations of DRAM output are generated 

DRAM v2 was implemented in [Nextflow](https://www.nextflow.io/) due to its innate scalability on HPCs and containerization, ensuring rigorous reproducibility and version control, thus making it ideally suited for high-performance computing environments. 

For more detail on DRAM and how DRAM v2 works please see our DRAM products
[DRAM version 1 publication](https://academic.oup.com/nar/article/48/16/8883/5884738)
[DRAM in KBase publication](https://pubmed.ncbi.nlm.nih.gov/36857575/)
[DRAM webinar](https://www.youtube.com/watch?v=-Ky2fz2vw2s)

### DRAM v2 Development Note

The DRAM development team is actively working on DRAM v2. We do not anticipate adding any additional functionality to DRAM v1.

## Quick links

- [Installation](#installation)
- [Requirements](#requirements)
- [General Instructions](#general-instructions)
- [Important Installation Notes](#important-computation-notes)
- [Important Computation Notes](#important-computation-notes)
- [Which is Better?](#which-is-better)
- [DRAM v2 Databases](#dram-v2-databases)
- [Example command-line usage](#example-usage)
- [All command-line options](#command-line-options)
- [Cool Nextflow Tips and Tricks](#nextflow-tips-and-tricks)
- [Resource Management](#resource-management)
- [Summary](#summary)
- [Citation](#citing-dram)

## Installation

### Requirements

* Nextflow >= v23.04.2.5870
* Some form of Conda or a Nextflow supported Container Runtime that can run Docker Images (Apptainer, Singularity CE, Docker, Podman, Sarus, etc.)
* Docker image file (if using a Container Runtime)
* DRAM databases (preformatted and downloaded via Globus, or with KEGG, formatted by the user)

### General Instructions:

1) If you do not have Nextflow installed, please follow the instructions [from their documentation page](https://www.nextflow.io/docs/stable/install.html).

2) Choose if you want to install DRAM v2 using Conda or a container runtime (Apptainer, Singularity CE, Docker, Podman, Sarus, etc.).
    - If you choose to use Conda, you will need to have Conda installed on your system. Miniconda is a good option for this because it is lightweight and easy to install. You can find the installation instructions [here](https://docs.anaconda.com/miniconda/miniconda-install/).
    - If you choose to use a container runtime, you will need to have the container runtime installed on your system. Apptainer is a good option for this, it is the modern and open-source version of Singularity. You can find the installation instructions [here](https://apptainer.org/docs/admin/main/installation.html) but any of the listed above should work.

3) Create a DRAM directory to store the DRAM files and the nextflow configuration file:

    ```
    mkdir DRAM
    cd DRAM
    ```

4) Download the DRAM data from Globus in this DRAM directory. This is a large download (> 500 GB) and will take some time. 
    - You can find the data [here](https://app.globus.org/file-manager?origin_id=97ed64b9-dea0-4eb5-a7e0-0b50ac94e889). UUID: 97ed64b9-dea0-4eb5-a7e0-0b50ac94e889 You will need a Globus account to access and download the data.
    - The download consists of a database folder that contains all the preformatted databases with the description database. As well as a containers folder that contains the docker image to be ran with your container runtime. 

5) DRAM by default looks for the databases and container image relative to the launch directory, if you would like to change this, change where you want to store the DRAM data, or other configuration options such as what container runtime you are using, SLURM options, customize or add other profile options, etc. you can download this defaults configuration file to customize your DRAM run.

    ```
    curl -o nextflow.config https://raw.githubusercontent.com/WrightonLabCSU/DRAM/refs/heads/dev/nextflow.config
    ```

    If you don't download the nextflow.config file, DRAM will run with the default settings, but you can still override options on the command line.

6) You can use nextflow to install and update DRAM by running the following command anywhere on your system:

    ```
    nextflow pull WrightonLabCSU/DRAM -r dev
    ```

    (Once DRAM v2 is out of development, the `-r dev` will be removed and the command will be `nextflow pull WrightonLabCSU/DRAM`)

    You can update DRAM by rerunning the pull command,

    ```
    nextflow pull DRAM
    ```

    You can specify a specific branch or version tag by adding `-r <branch or tag>` to the pull command. For example, to pull the development branch, you would run (Note only version tags >= 2 are supported, pre 2 versions are before nextflow was used):

    ```
    nextflow pull WrightonLabCSU/DRAM -r dev
    ```

    Then to run DRAM, you can use the following command:

    ```
    nextflow run DRAM <options>
    ```
   
    You can update DRAM by rerunning the pull command,


### Important Installation Notes:

Nextflow installs all nextflow pipeline scripts by default in the `$HOME/.nextflow/assets/WrightonLabCSU/DRAM` directory, allowing nextflow to centralize install and update management. If you are running DRAM on a shared system, you may wish to install DRAM in a shared directory. DRAM's actual pipeline scripts are very small and if multiple users are running DRAM on the same system, having them each install their own copy of DRAM does not take up a lot of space, but it does make it easier to manage updates and versions. 

Two ways to install DRAM in a shared directory are:

1) Where nextflow installs all of its scripts when you do `nextflow pull` is defined in the environment variable `NXF_HOME`, which if not set defaults to `$HOME/.nextflow`. You can set this on your system to a shared directory such as `/home/nextflow` and then run `nextflow pull` to install DRAM in that shared directory. This will install all nextflow scripts in that shared directory. You may or may not want to install all nextflow scripts in a shared directory, so this may not be the best option.

2) In the launch directory `nextflow.config`, you can add a section that tells nextflow to set `NXF_HOME` for the DRAM pipeline. This will only set `NXF_HOME` for the DRAM pipeline and not for all nextflow pipelines. This is the best option if you only want to set `NXF_HOME` for the DRAM pipeline. At the top of your `nextflow.config` file, you could add something like:

    ```
    env {
        NXF_HOME='$PWD'
    }
    ```

    This will set `NXF_HOME` to the CWD of the config file for the DRAM pipeline only (or set it anywhere else you want).


You will be able to run DRAM from anywhere, not just the directory you downloaded your data and config too with the command `nextflow run DRAM`, but if it is anywhere other than where the nextflow.config file is located, you will need to specify the path to the nextflow.config file with the `-c` flag.

### Important Computation Notes:

DRAM utilizes either Conda or Singularity for dependency management and the user MUST choose one of the following options on execution of any DRAM command

*The Nextflow profile option is used (`-profile`) - yes! a single hyphen! Nextflow options use a single hyphen, while DRAM options use the traditional double hyphen*

1) `-profile conda`
   
  This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

2) `-profile conda_slurm`
   
  This option will submit each individual DRAM process as its own SLURM job. (See Wiki Resource Management for details).
  This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

3) `-profile singularity`
   
  This option relies on the local systems Singularity. Nextflow will create its own Conda environments to run in. 

4) `-profile singularity_slurm`
   
  This option will submit each individual DRAM process as its own SLURM job. (See Wiki Resource Management for details)
  This option relies on the local systems Singularity to run the downloaded Singularity container.  


### Which is Better?

#### Conda Environments

#### Pros:

- Beginner-Friendly: Easy to install and use, making it accessible for newcomers.
- Reproducibility: Efficient management of environments facilitates reproducibility.

##### Cons:

- Generally slower than using Singularity containers (will have metrics in the future).
- Dependency Conflicts: Dependency resolution can be slow and may lead to conflicts.
- Limited Portability: System dependencies may introduce variability, affecting portability.
- System Variability: Reliance on the host system's architecture and libraries can cause variability between systems.

#### Singularity Containers

##### Pros:

- Generally faster than using Conda environments (will have metrics in the future).
- Consistent Environments: Ensures consistent runtime environments, enhancing reproducibility.
- HPC Ideal: Perfect for high-performance computing (HPC) environments without the need for root access.
- Isolation: Offers isolation from the host system, minimizing conflicts.
- Wide Portability: Containers are portable across any Linux system with Singularity.

##### Cons:

- Installation Complexity: Can be trickier to install compared to Conda environments.
- Storage Space: May consume more storage space.

Conda is recommended for its ease of use and versatility across different programming languages.
Singularity excels in ensuring reproducibility and compatibility in high-performance computing environments.

---------

## DRAM v2 Databases

DRAM v2 databases, unlike DRAM v1 databases, will be pre-formatted and hosted online. Users of DRAM will need to:
1) decide which databases suits their needs
2) download the databases from Globus

  *However, these databases can be quite large and it is therefore important to look through the options below.*

  *These databases rely on an SQL database of database descriptions which is provided in 3 different sizes based on ther user's needs.*

#### All databases DRAM accommodates:
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

<a name="database-sets"></a>


## Example usage

DRAM apps Call, Annotate and Distill can all be run at once or alternatively, each app can be run individually (assuming you provide the required input data for each app).

Additionally, `--merge-annotations` and `--rename` can be run idenpendently of any other apps. 


1) **Rename fasta headers based on input sample file names:**

`nextflow run DRAM --rename --input_fasta_dir <path/to/fasta/directory/>`


2) **Call genes using input fastas (use --rename to rename FASTA headers):**

`nextflow run DRAM --call --rename --input_fasta_dir <path/to/fasta/directory/>`

    
3) **Annotate called genes using input called genes and the KOFAM database:**

`nextflow run DRAM --annotate --input_genes <path/to/called/genes/directory> --use_kofam`


4) **Annotate called genes using input fasta files and the KOFAM database:**

`nextflow run DRAM --annotate --input_fasta <path/to/called/genes/directory> --use_kofam`


5) **Merge verious existing annotations files together (Must be generated using DRAM.):**

`nextflow run DRAM --merge_annotations <path/to/directory/with/multiple/annotation/TSV/files>`


6) **Distill using input annotations:**

`nextflow run DRAM --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>`


7) **(Combined): Call, annotate and distill input fasta files:**

`nextflow run DRAM --rename --call --annotate --use_<database(s) --distill_topic <distillate(s)>`


8) **Call and Annotate genes using input fastas and KOFAM database. Distill using the default topic and the AG ecosystem:**

`nextflow run DRAM --input_fasta_dir <path/to/fasta/directory/> --outdir <path/to/output/directory/> --call --annotate --distill_topic default --distill_ecosystem ag --threads <threads> --use_kofam`


9) **"Real-world" example using the test data provided in this repository:**
   
(put on multiple lines for clarity)

```
nextflow run -bg
DRAM
--input_fasta ../test_data/DRAM_test_data/
--outdir DRAM-test-data-call-annotate-distill
--threads 8
--call --rename --annotate
--use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper
--use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur
--distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom test-data/custom-test-distilalte.tsv
-profile conda_slurm --slurm_node main -with-report -with-trace -with-timeline
```

  **Breakdown of example (9):**
  - `-bg` : Nextflow option to push the run immediately into the background. (Thus, you can log out on an HPC and the run will continue).
  - `-profile` : Nextflow option to select profile (Conda vs Singularity and SLURM vs no-SLURM).
  - `--slurm_node`: DRAM option to select a specific node to compute on during the whole run.
  - `-with-trace`: Nextflow option to output a process-by-process report of the run. (TEXT)
  - `-with-report`: Nextflow option to output a process-by-process report of the run. (HTML)
  - `-with-timeline`: Nextflow option to output a process-by-process HTML timeline report of the run. (HTML)

-------

## Command-line Options

### General Command-line Options

    Description: 
        The purpose of DRAM is to provide FASTA annotation, across a vast array of databases, with expertly-currated distillation. 
        DRAM can be used to call, annotate and distill annotations from input FASTA files. 
        Call, annotate and distill can be run together or, each can be run idependently. 

    Bring up help menu:
        nextflow run DRAM --help (--h)

    Bring up versions menu:
        nextflow run DRAM --version (--v)      

    Usage:
        nextflow run DRAM --rename --call --annotate --use_<database(s) --distill_topic <distillate(s)>

        Call genes using input fastas (use --rename to rename FASTA headers):
            nextflow run DRAM --call --rename --input_fasta_dir <path/to/fasta/directory/>

        Annotate called genes using input fastas:
            nextflow run DRAM --annotate --input_genes <path/to/called/genes/directory>

        Distill using input annotations:
            nextflow run DRAM --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>

        (Combined): Call, annotate and distill input fasta files:
            nextflow run DRAM --rename --call --annotate --use_<database(s) --distill_topic <distillate(s)>

        (Real) example: (on multiple lines for clarity)
        nextflow run DRAM --input_fasta ../test_data/ 
            --outdir DRAM-test-data-Feb012024/ 
            --call --rename 
            --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur 
            --add_annotations ../test-data/old-DRAM-annotations.tsv
            --distill_topic 'carbon transport energy' --distill_ecosystem 'eng_sys ag' 
            --distill_custom assets/forms/distill_sheets/test.tsv -resume --slurm_node zenith 
            --trnas ../test-data/trnas.tsv
            --rrnas ../test-data/rrnas.tsv
            --bin_quality ../test-data/checkM1-test-data.tsv
            --taxa ../test-data/gtdbtk.bac120.summary.tsv
            --generate_gff 
            --generate_gbk
            --threads 5
            -with-report -with-trace -with-timeline

    Main DRAM Operations:
        --call      : Call genes using prodigal 
        --annotate  : Annotate called genes using downloaded databases
        --distill   : Distill the annotations into a multi-sheet distillate.xlsx
    
    REQUIRED DRAM profile options:
        -profile                STRING  <conda, conda_slurm, singularity, singularity_conda>
                                        Runs DRAM either using Conda (must be installed) or Singularity (must be installed).
                                        Runs DRAM with no scheduling or scheduling via SLURM.
                                        See SLURM options in full help menu.

    Call options:
        --call                  OPTION  Call genes on the input FASTA files using Prodigal.

        --input_fasta           PATH    <path/to/fasta/directory/>
                                            Directory containing input fasta files.      
                                            Default: <./input_fasta/*.fa*>

        --rename                OPTION  Rename FASTA headers based on file name.    
                                            Example: sample1.fa --> (fasta header renamed to) > sample1......
                                            Why? DRAM output is focused on scaffolds/contigs with respect to each provided input sample.
                                                Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguashable.
                                                *If you have already renamed your FASTA headers, do not include '--call'.

        --prodigal_mode         STRING  <single|meta>
                                            Default: 'single'

        --prodigal_tras_table   NUMBER  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
                                            Specify a translation table to use (default: '1').
                                    
        --min_contig_len        NUMBER  <number in base pairs>
                                            Default: '2500'
                                            
    Annotate options:
        --use_<db-name>         STRING   <camper|cant_hyd|dbcan|fegenie|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref]
                                            Specify databases to use. Can use more than one. Can be used in combination with --use_dbset.
        
        --use_dbset             STRING  <metabolism_kegg_set|metabolism_set|adjectives_kegg_set|adjectivs set>
                                            metabolism_kegg_set = kegg, dbcan, merops, pfam, heme
                                            metabolism_set      = kofam, dbcan, merops, pfam, heme
                                            adjectives_kegg_set = kegg, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                            adjectives_set      = kofam, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                            *Only one set can be used. Can be used in combination with --use_[db-name]
        
        --input_genes           PATH    <path/to/called/genes/directory/>
                                            Directory containing called genes (.faa) 

        --add_annotations       PATH    <path/to/old-annoations.tsv> 

        --generate_gff          OPTION Will generate an output GFF for each sample based on the raw-annotations.tsv.

        --generate_gbk          OPTION Will generate an output GBK for each sample based on the raw-annotations.tsv.
                                                    Used to add in old annotations to the current run. (See example for format.)

    Distill options:
        --annotations           PATH     <path/to/annotations.tsv>
                                            Required if you are running distill without --call and --annotate.

        --rrnas                 PATH    <path/to/rRNA.tsv> (See example for format.)
                                            rRNA information will be included in distill output.

        --trnas                 PATH    <path/to/tRNA.tsv> (See example for format.)
                                            tRNA information will be included in distill output.
        
        --bin_quality           PATH    <path/to/bin-quality.tsv> (See example for format.)
                                            CheckM and CheckM2 compatible. 

        --taxa                  PATH    <path/to/bin-taxonomy.tsv>
                                            Compatible with GTDB. (See example for format.)

        --distill_topic         STRING  <carbon|energy|misc|nitrogen|transport> OR <default = carbon, energy, misc, nitrogen, transport>
                                            If more than one topic included, they must be enclosed in single quotes

        --distill_ecosystem     STRING  <eng_sys|ag>
                                            If more than one ecosystem included, they must be enclosed in single quotes

        --distill_custom        STRING  <path/to/custom_distillate.tsv> (See example for format and options.)
                                            As of now, only one custom distillate may be included.

    General options:
        --outdir                PATH    <path/to/output/directory>
                                            Default: './DRAM_output/'

        --threads               NUMBER  Number of threads to use for processing.
                                        Default: '10'

        --slurm_node            STRING  <node_name>
                                        Example --slurm_queue c001

        --slurm_queue           STRING  <slurm partition name>
                                        Example:  --slurn_queue 'smith-hi,smith-low'

        -with-trace             OPTION  Nextflow option to output a process-by-process report of the run. (TEXT)
        -with-report            OPTION  Nextflow option to output a process-by-process report of the run. (HTML)
        -with-timeline:         OPTION  Nextflow option to output a process-by-process HTML timeline report of the run. (HTML)

### Call Command-line Options

    Call description: The purpose of DRAM --call is to call genes on input FASTA files.

    Usage:

        Call genes using input fastas:
            nextflow run DRAM --call --input_fasta_dir <path/to/fasta/directory/> --outdir <path/to/output/directory/> --threads <threads>
    
    REQUIRED DRAM profile options:
        -profile                STRING  <conda, conda_slurm, singularity, singularity_conda>
                                        Runs DRAM either using Conda (must be installed) or Singularity (must be installed).
                                        Runs DRAM with no scheduling or scheduling via SLURM.
                                        See SLURM options in full help menu.

    Call options:
        --rename                Rename FASTA headers based on file name.    
                                    Example: sample1.fa --> (fasta header renamed to) > sample1......
                                    Why? DRAM output is focused on scaffolds/contigs with respect to each provided input sample.
                                        Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguashable.
                                        *If you have already renamed your FASTA headers, do not include '--call'.

        --prodigal_mode         STRING  <single|meta>
                                    Default: 'single'

        --prodigal_tras_table   <1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25>
                                    Specify a translation table to use (default: '1').
                                    
        --min_contig_len        NUMBER  <number in base pairs>
                                            Default: '2500'

    Main options:
        --input_fasta           PATH    <path/to/fasta/directory/>
                                        Directory containing input fasta files.      
                                        Default: './input_fasta/'

        --outdir                PATH    <path/to/output/directory>
                                            Default: './DRAM_output/'

        --threads               NUMBER  Number of threads to use for processing.
                                        Default: '10'

        --slurm_node            STRING  <node_name>
                                        Example --slurm_queue c001

        --slurm_queue           STRING  <slurm partition name>
                                        Example:  --slurn_queue 'smith-hi,smith-low'

        -with-trace             OPTION  Nextflow option to output a process-by-process report of the run. (TEXT)
        -with-report            OPTION  Nextflow option to output a process-by-process report of the run. (HTML)
        -with-timeline:         OPTION  Nextflow option to output a process-by-process HTML timeline report of the run. (HTML)
                                       
### Annotate Command-line Options

    Annotate description: The purpose of DRAM '--annotate' is to annotate called genes on input (nucleotide) FASTA (fa*) files.

    Usage:

        Annotate called genes using input called genes and the KOFAM database:
            nextflow run DRAM --annotate --input_genes <path/to/called/genes/directory> --use_kofam
        
        Annotate called genes using input fasta files and the KOFAM database:
            nextflow run DRAM --annotate --input_fasta <path/to/called/genes/directory> --use_kofam
    
    REQUIRED DRAM profile options:
        -profile            STRING  <conda, conda_slurm, singularity, singularity_conda>
                                        Runs DRAM either using Conda (must be installed) or Singularity (must be installed).
                                        Runs DRAM with no scheduling or scheduling via SLURM.
                                        See SLURM options in full help menu.

    Annotate options:
    --use_<db-name>         STRING   <camper|cant_hyd|dbcan|fegenie|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref]
                                        Specify databases to use. Can use more than one. Can be used in combination with --use_dbset.
    
    --use_dbset             STRING  <metabolism_kegg_set|metabolism_set|adjectives_kegg_set|adjectivs set>
                                        metabolism_kegg_set = kegg, dbcan, merops, pfam, heme
                                        metabolism_set      = kofam, dbcan, merops, pfam, heme
                                        adjectives_kegg_set = kegg, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                        adjectives_set      = kofam, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
                                        *Only one set can be used. Can be used in combination with --use_[db-name]
    
    --add_annotations       PATH    <path/to/old-annoations.tsv> 
                                        Used to add in old annotations to the current run. (See example for format.)

    --generate_gff          OPTION Will generate an output GFF for each sample based on the raw-annotations.tsv.

    --generate_gbk          OPTION Will generate an output GBK for each sample based on the raw-annotations.tsv.
    
    Main options:
    --input_fasta           PATH    <path/to/fasta/directory/>
                                        Directory containing input fasta files.      
                                        Default: './input_fasta/' 
                                        Either '--input_fasta' or '--input_genes' may be used - not both.

    --input_genes           PATH    <path/to/called/genes/directory/>
                                        Directory containing called genes (.fna)
                                        Either '--input_fasta' or '--input_genes' may be used - not both.

    --outdir                PATH    <path/to/output/directory/>
                                        Default: './DRAM_output/'

    --threads               NUMBER  Number of threads to use for processing.
                                        Default '10'

    --slurm_node            string  <node_name>
                                    Example --slurm_queue c001

    --slurm_queue           string  <slurm partition name>
                                    Example:  --slurn_queue 'smith-hi,smith-low'

    -with-trace             OPTION  Nextflow option to output a process-by-process report of the run. (TEXT)
    -with-report            OPTION  Nextflow option to output a process-by-process report of the run. (HTML)
    -with-timeline:         OPTION  Nextflow option to output a process-by-process HTML timeline report of the run. (HTML)         
                                
### Distill Command-line Options

    DRAM Nextflow Pipeline
    ===================================
    Distill description:    The purpose of DRAM --distill is to distill down annotations based on curated distillation summary form(s). 
                            User's may also provide a custom distillate via --distill_custom <path/to/file> (TSV forms).
                            Distill can be ran independent of --call and --annotate however, annotations must be provided (--annotations <path/to/annotations.tsv>). 
                            Optional tRNA, rRNA and bin quality may also be provided.
    
    Usage:
        nextflow run DRAM --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/> --threads <threads>
        *Important: if more than one topic or ecosystem is included, they must be enclosed in single quotes. Example: --distill_topic 'carbon transport'
    
    Example:
        Call and Annotate genes using input fastas and KOFAM database. Distill using carbon topic and AG ecosystem:
            nextflow run DRAM --input_fasta_dir <path/to/fasta/directory/> --outdir <path/to/output/directory/> --call --annotate --distill_topic carbon --distill_ecosystem ag --threads <threads> --use_kofam
    
    REQUIRED DRAM profile options:
        -profile                STRING  <conda, conda_slurm, singularity, singularity_conda>
                                        Runs DRAM either using Conda (must be installed) or Singularity (must be installed).
                                        Runs DRAM with no scheduling or scheduling via SLURM.
                                        See SLURM options in full help menu.

    Distill options:
        --annotations           PATH     <path/to/annotations.tsv>
                                            Required if you are running distill without --call and --annotate.

        --rrnas                 PATH    <path/to/rRNA.tsv> (See example for format.)
                                            rRNA information will be included in distill output.

        --trnas                 PATH    <path/to/tRNA.tsv> (See example for format.)
                                            tRNA information will be included in distill output.

        --bin_quality           PATH    <path/to/bin-quality.tsv> (See example for format.)
                                            CheckM and CheckM2 compatible. 

        --taxa                  PATH    <path/to/bin-taxonomy.tsv>
                                        Compatible with GTDB. (See example for format.)

        --distill_topic         STRING  <carbon|energy|misc|nitrogen|transport> OR <default = carbon, energy, misc, nitrogen, transport>
                                            If more than one topic included, they must be enclosed in single quotes

        --distill_ecosystem     STRING  <eng_sys|ag>
                                            If more than one ecosystem included, they must be enclosed in single quotes

        --distill_custom        STRING  <path/to/custom_distillate.tsv> (See example for format and options.)
                                            As of now, only one custom distillate may be included.

    Main options:                
        --outdir                PATH    <path/to/output/directory/>
                                            Default: './DRAM_output/'

        --threads               NUMBER  Number of threads to use for processing.
                                            Default '10'

        --slurm_node            string  <node_name>
                                        Example --slurm_queue c001

        --slurm_queue           string  <slurm partition name>
                                        Example:  --slurn_queue 'smith-hi,smith-low'

        -with-trace             OPTION  Nextflow option to output a process-by-process report of the run. (TEXT)
        -with-report            OPTION  Nextflow option to output a process-by-process report of the run. (HTML)
        -with-timeline:         OPTION  Nextflow option to output a process-by-process HTML timeline report of the run. (HTML)
-----------
------

### Nextflow Tips and Tricks

In Nextflow DSL2, the `-resume` option offers a powerful feature that allows you to efficiently manage and modify your workflow runs. It enables you to resume a run after it has finished, make changes to parameters, and reuse previously generated data, enhancing the flexibility and reusability of your workflow. Here are some common scenarios where the `-resume` option comes in handy:

- **I Want to run DRAM again but with additional databases**
  - By using the `-resume` option, if you have not deleted your `work/` directory, running your previous command with the desired changes will allow you to reuse your called genes and existing annotations.
  - For example, if you intially used `--use_kofam --use_dbcan`, you can add in additional annotation databases like, `--use_kegg --use_uniref`. This will reuse the existing called genes to annotate the newly added databases and will retain the existing annotations (assuming they were retained within the modified command).

These Nextflow tips and tricks demonstrate how the `-resume` option can optimize your workflow, save time, and improve the reusability of previously computed data.

**How the resume option does NOT work:**
- <examples> TODO: Add examples of how it won't work

Every CLI option in DRAM shown above or from `nextflow run DRAM --help` can be set in the `nextflow.config` file. This can be useful for setting default options or for setting options that are the same for every run. For example, if you always want to use the `--use_uniref` option, you can set this in the `nextflow.config` file and it will be used for every run by setting in the params section of the `nextflow.config` file:

    ```
    params {
        use_uniref = true
        // or use_uniref = 1
    }
    ```

If you are always running `--annotate`, you can set this in the `nextflow.config` file and it will be used for every run by setting in the params section of the `nextflow.config` file:

    ```
    params {
        annotate = true
        // or annotate = 1
    }
    ```


Nextflow will use the `nextflow.config` file in the launch directory by default. If you want to use a different `nextflow.config` file, you can specify it on the command line with the `-c` option:

    ```
    nextflow run DRAM -c /path/to/some_config_file.config
    ```

If this is done, the `some_config_file.config` file will only override the options in the `nextflow.config` file in the launch directory that are defined in the `some_config_file.config` file. All other options will be taken from the `nextflow.config` file in the launch directory. 

This can be useful if you have a set of options that you want to use for a specific run, but you want to keep the default options for all other runs. This can also be useful in a shared user environment where you want to keep your options separate from other users.

--------

### Resource Management

DRAM is implemented in Nextflow.

Nextflow is able to scale horizontally on a given machine.

What does this mean?

Horizontal scaling refers to the ability to distribute computational tasks across multiple computing resources, such as cores or nodes, in parallel. In the context of Nextflow, this means that a single workflow can leverage the computational power of multiple CPUs or nodes, allowing for faster execution of tasks and improved overall performance.

By utilizing horizontal scaling, Nextflow can efficiently manage and execute workflows that require significant computational resources, such as those involved in genomic data analysis. This enables DRAM to process large datasets and complex analyses in a timely manner, making it suitable for a wide range of research and bioinformatics applications.

## Summary

DRAM comes with configuration files which have the option to change how many "things" can happen at a time in the pipeline.

A user can modify these, "maxForks", parameters within the ``nextflow.config`` to increase the number of "things" which DRAM can perform at a given time.

**NOTE**: Development is in progress to enable different DRAM modes: "lite", "medium" and "heavy". Where "lite" would be for a good laptop and "heavy" for a HPC. These options will alter the CPU and memory (RAM) requirements for each process.

---------------

## Citing DRAM
The DRAM was published in Nucleic Acids Research in 2020 and is available [here](https://academic.oup.com/nar/article/48/16/8883/5884738). If
DRAM helps you out in your research, please cite it.
