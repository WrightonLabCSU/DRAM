<p align="center">
  <img src="assets/images/DRAM2_large.png" width="600" height="600" alt="DRAM v2 logo">
</p>

## ⚠️ DRAM v2 is currently under active development and usage is at your own risk. ⚠️


DRAM v2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic and genomic assembled data (e.g. scaffolds or contigs) or called genes (e.g. nuclotide or amino acid format). DRAM annotates MAGs using [KEGG](https://www.kegg.jp/) (if provided by the user), [UniRef90](https://www.uniprot.org/), [PFAM](https://pfam.xfam.org/), [dbCAN](http://bcb.unl.edu/dbCAN2/), [RefSeq viral](https://www.ncbi.nlm.nih.gov/genome/viruses/), [VOGDB](http://vogdb.org/) and the [MEROPS](https://www.ebi.ac.uk/merops/) peptidase database as well as custom user databases. DRAM is run in two stages. First an annotation step to assign database identifiers to gene, and then a distill step to curate these annotations into useful functional categories. DRAM was implemented in [Nextflow](https://www.nextflow.io/) due to its innate scalability on HPCs and containerization, ensuring rigorous reproducibility and version control, thus making it ideally suited for high-performance computing environments. 

DRAM is run in four stages: 
1) Gene Calling Prodogal - genes are called on user provided scaffolds or contigs 
2) Gene Annotation - genes are annotated with a set of user defined databases 
3) Distillation - annotations are curated into functional categories
4) Product Generation - interactive visualizations of DRAM output are generated 

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
- [DRAM v2 Databases](#dram-v2-databases)
- [Example command-line usage](#example-usage)
- [All command-line options](#command-line-options)
- [Cool Nextflow Tips and Tricks](#nextflow-tips-and-tricks)
- [Resource Management](#resource-management)
- [Summary](#summary)
- [Citation](#citing-dram)

## Installation

### Requirements

* Nextflow >= v23.04.2
* Some form of Conda or a Nextflow supported Container Runtime (Apptainer, Singularity CE, Docker, Podman, Sarus, etc.)
* DRAM databases (preformatted and downloaded via Globus, or with KEGG, formatted by the user)

### General Instructions:

On many HPC systems, Nextflow, Conda, and some container runtime (usually Singularity or Apptainer) are already installed. If you are using a local system, you will need to install Nextflow and Conda or a container runtime. On an HPC system, Nextflow, Conda, and a container runtime are often loaded with modules. Refer to your HPC system's documentation for instructions on how to load modules. They might be loaded another way or already pre-installed.

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
    - The download consists of a database folder that contains all the preformatted databases with the description database.

5) DRAM by default looks for the databases relative to the launch directory, if you would like to change this, change where you want to store the DRAM data, or other configuration options such as what container runtime you are using, SLURM options, customize or add other profile options, etc. you can download this defaults configuration file to customize your DRAM run.

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

The best way to accomplish this is to download or clone the DRAM repository to a shared directory and launch the `main.nf` script in the root directory with the following command:

```
nextflow run <path/to/DRAM>/main.nf <options>
```

You can use the `-c` to specify the path to a custom nextflow.config file. For example, if you have a custom nextflow.config file in another directory, you can run the following command:

```
nextflow run <path/to/DRAM>/main.nf -c <path/to/custom/nextflow.config> <options>
```

### Important Computation Notes:

DRAM comes with a variety of profiles to choose from. The profiles are used to specify the environment for dependency management. The main options are conda, singularity, apptainer, mamba, docker, podman, and shifter. Others can be found in the `nextflow.config` file in the root of the DRAM repo. It is specified with `-profile <profile_name>`.

*The Nextflow profile option is used (`-profile`) - yes! a single hyphen! Nextflow options use a single hyphen, while DRAM options use the traditional double hyphen*

Use the `-profile` option that best fits your needs. Many HPC systems come with a way to load modules that give access to conda, apptainer, and singularity. The initial time you load a profile, the dependencies will be installed and then cached in a temp directory called `work` in the launch directory (See the nextflow docs page for information on specifing the cache directory directly). Then it will be reused for subsequent runs.

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

Additionally, `--merge_annotations` and `--rename` can be run idenpendently of any other apps. 

**You can also run DRAM with the `--slurm` option to have nextflow launch jobs on a SLURM cluster. The actual nextflow job is lightweight enough to usually run as a background process in your main session. If you are launching nextflow jobs from you main session on a server where you might get disconnected, launch you run with the `-bg` flag to run it as a background task so it doesn't get killed when you get disconnected**


1) **Rename fasta headers based on input sample file names:**

`nextflow run DRAM --rename --input_fasta <path/to/fasta/directory/>`


2) **Call genes using input fastas (use --rename to rename FASTA headers):**

`nextflow run DRAM --call --rename --input_fasta <path/to/fasta/directory/>`

    
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

`nextflow run DRAM --input_fasta <path/to/fasta/directory/> --outdir <path/to/output/directory/> --call --annotate --distill_topic default --distill_ecosystem ag --threads <threads> --use_kofam`


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

See `--help` for a full list of options. A docs page with all the options is in progress.

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
    }
    ```

If you are always running `--annotate`, you can set this in the `nextflow.config` file and it will be used for every run by setting in the params section of the `nextflow.config` file:

    ```
    params {
        annotate = true
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
