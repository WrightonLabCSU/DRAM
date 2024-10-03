# DRAM2 Singularity containers and Conda environments

Here will be described how and why DRAM2 relies on Singularity containers and, ideally as a backup, Conda environments. At the end are a few future development notes.


## Why Singularity containers

Containerization is a technology that allows software applications and their dependencies to be packaged together into a single, portable unit called a container. Containers ensure consistency and stability by encapsulating everything the application needs to run, making it easier to deploy across various environments without compatibility issues. In Nextflow, containers are often used to run pipeline processes, ensuring that each process has the exact environment it requires, which enhances reproducibility and version control. By using containerization platforms like Docker or Singularity, Nextflow workflows can be executed reliably on different systems, from local machines to cloud-based infrastructures. This approach minimizes the risk of software conflicts and simplifies the maintenance of complex bioinformatics workflows.

------

## Why Conda containers (ideally as a backup to Singularity containers)

Conda environments can also be used in Nextflow pipelines to manage dependencies and ensure that software packages are correctly installed. Conda is particularly useful because it allows for easy creation and management of isolated environments with specific package versions. However, it is recommended to use Conda as a backup to Singularity containers because Conda can sometimes be tricky to manage across different system architectures. Issues such as inconsistencies in package availability and compatibility can arise when using Conda on diverse systems. Therefore, while Conda is a flexible tool, Singularity containers provide a more stable and reproducible solution for running Nextflow workflows.

------

## What are Singularity containers, how are they built and what containers are used for DRAM2

### What are Singularity containers

Singularity containers are essentially a mini operating system in which dependencies are installed. A user writes commands into a recipe file, just like a user would write commands into the command-line to install dependencies. Then, a container is built based on the recipe file. If the container is built successfully, this container is now a self-contained environment that includes all the necessary software and libraries. This container can be moved across different systems, ensuring that the software will run exactly the same way regardless of the underlying hardware or operating system. Singularity is particularly well-suited for high-performance computing environments because it maintains a high level of security and integration with the host system, making it an ideal choice for reproducible research and complex workflows in Nextflow.

### How are Singularity containers built and how are they referenced in DRAM2

To build a Singularity container, you can use the following command, where `DRAM2-Nextflow-Main-Container-June282024-V5` is the name of your recipe file:

```bash
sudo singularity build DRAM2-Nextflow-Main-Container-June282024-V5.sif DRAM2-Nextflow-Main-Container-June282024-V5
```

*Note: `sudo` is used here to build the container. While Singularity does have a `--fake-root` option, it is advised to use `sudo`. If a user does not have `sudo` privileges on a machine, it is advised to build them locally, where the user has sudo, and upload them to the machine used for DRAM2 development.*

This command will create a .sif file, which is the Singularity Image Format, containing all the dependencies and software specified in your recipe file. This image can then be used to run processes in your Nextflow pipeline, ensuring a consistent and reproducible environment.

In the `nextflow.config` file on line 296 we define `main_container` as the built `.sif` container:

```bash
    /* Containers and Environments */
        main_container = "./containers/DRAM2-Nextflow-Main-Container-June282024-V5.sif"
```

Then, `main_container` is subsequently used within the profile definitions:

Example of `singularity_slurm` profile referencing the `main_container`:
```bash
profiles {
    <...>
    singularity_slurm {
        includeConfig "${params.singularity_slurm_config}"
        process.container = "${params.main_container}"
        singularity {
            enabled = true
            autoMounts = true
        }
        executor{
            queueSize = "${params.queue_size}"
        }
    }
```

### DRAM2 containers

DRAM2 only has one container! 

DRAM2 `main_container`:
`DRAM2-Nextflow-Main-Container-March262024-V5.sif`

#### Why only one container

This is ideal as images can be large in size. As a software grows, and dependencies become numerous, it can be difficult for all of the dependencies to reside within a single container. However, this can usually be solved by creating Conda environments within a given container. Additionally, installing Conda within a container is typically an easier way to install dependencies within a given container.

Here is an example of installing Conda within a container and performing subsequent installs:

```bash
    # Download and install Miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda
    export PATH="/opt/miniconda/bin:$PATH"

    # Add conda channels
    conda config --add channels anaconda
    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels agbiome
    conda config --add channels r
    conda update -n base -c defaults conda

    conda install -c bioconda prodigal=2.6.3
```

-------

## Conda environments

As stated above, Conda environment are ideally used as a backup. However, their use, construction, and implementation will be described here. 

While only one Singularity container is used, multiple Conda environments are needed in the event the user is unable to use Singularity. 

### How are Conda environments built and how are they referenced in DRAM2

#### Building Conda environments - writing recipe files

Conda environments do NOT have to be built before running DRAM2 by a given user. 

**This does not mean they do not have to be tested.**

Nextflow will use a provided Conda recipe file `*.yml` file to build Conda containers as they are needed upon pipeline execution. Additionally, if a user runs multiple DRAM2 runs, and they retain the `work/` directory, Nextflow will rely on cached built Conda environments.

Thus, recipe files must be provided with the required dependencies. 

Here is an example of the `environment-quast.yml` Conda recipe file:

```bash
name: dram2-quast-env
channels:
  - conda-forge
  - bioconda
  - defaults
  - agbiome

dependencies:
  - python>=3.8
  - pandas
  - openpyxl=3.1.2
  - quast=5.2.0
  - pip:
      - matplotlib==3.8.3
```

*Note: QUAST for some reason is difficult to install within the `main-environment.yml` and thus there are two `*.yml` recipe files.

#### How are Conda environment recipe files referenced in DRAM2

Similar to Singularity containers, Conda recipe files are referenced in the `nextflow.config`

Example `nextflow.config` at line 297:

```bash
    /* Containers and Environments */
        main_container = "./containers/DRAM2-Nextflow-Main-Container-March252024-V4.sif"
        main_environment = "./assets/conda/environment.yml"
        quast_environment = "./assets/conda/environment-quast.yml"
```

Here are, the variables `main_environment` and `quast_environment` are defined.

Then, these vairiable are referenced within the definition.

Example Conda profile definition:

```bash
profiles {
    conda_slurm {
        includeConfig "${params.conda_slurm_config}"
        conda.enabled = true
        process.conda = "${params.main_environment}"
        executor{
            queueSize = "${params.queue_size}"
        }
    }
```
**NOTE: `process.conda` needs to be defined as a given container BUT this is purely to satisfy Nextflow profile syntax. One might be confused as why the `quast_environment` is not referenced within the `conda_slurm` profile definition. This is explained below.**

#### Referencing individual Conda environment within the individual profile configuration files.

As described in `DRAM2-Nextflow-Profiles.md` and parially, in `DRAM2-Computational-Resource-Management.md`, profiles are used to run DRAM2 using Singularity or Conda, with or without SLURM job submission.

Within a given profile configuration file, each DRAM2 process is defined and the computational resources are defined. 

Example portion of profile for Conda + SLURM:

```bash
process {
    withName: CALL_GENES {
        publishDir = [
            path: "${params.outdir}/Prodigal_v2.6.3",
            mode: 'copy',
        ]

        conda = "${params.main_environment}"
        executor="slurm"
        cpus = "${params.single_thread_job}"
        time = "${params.time}"
        memory = "${params.sml_job}"
        queue = "${params.slurm_queue}"
        clusterOptions = (params.slurm_node ? "--nodelist=${params.slurm_node} " : "") + "--job-name=DRAM2-call"

        maxForks = "${params.max_forks_single_cpu}"
    }

    withName: QUAST {
        publishDir = [
            path: "${params.outdir}/QUAST_v5.2.0",
            mode: 'copy',
        ]

        conda = "${params.quast_environment}"
        executor="slurm"
        cpus = "${params.threads}"
        time = "${params.time}"
        memory = "${params.sml_job}"
        queue = "${params.slurm_queue}"
        clusterOptions = (params.slurm_node ? "--nodelist=${params.slurm_node} " : "") + "--job-name=DRAM2-quast"

        maxForks = "${params.max_forks_single_cpu}"
    }
```

Here we have the processes CALL_GENES() and QUAST(), each of which use a different Conda recipe file.

QUAST uses:
```bash
conda = "${params.quast_environment}"
```

While CALL_GENES() uses:
```bash
conda = "${params.main_environment}"
```

--------

## Future development notes

From development of COMET, it was noticed that relying on Nextflow to build the Singularity containers resulted in conflicts across different HPCs, that is W2 and Riviera. Thus, it was decided, in COMET, to pre-build the containers. Within Nextflow you can, instead of linking the `*.yml` recipe files, provide the location to a pre-build Conda environment. However, this means either the user must build these (As is instructed for COMET) or DRAM2 must provide these.

It is suggested to store the Singularity containers, and pre-built Conda environments if that is chosen, on GLOBUS with the databases. Then a given user can pull everything at once. 

For future building of Singularity environments, for speed of building, it is suggested to install Mamba within the recipe files and use Mamba for package installation.

Once DRAM2 is released, and even before, it needs to be ensure that the dependency versions between the Singularity containers and the Conda recipe files match. Or, if this is not possible due to Conda difficulties, the differences must be stated in the help menus and in all documentation.

Lastly, Singularity is not as widespread as Docker. It is suggested to create Docker-based profiles and user the built-in Singularity functionality to convert the Singularity containers to Docker containers.


## General Singularity Container information for W2 computational team

**Note: This is general information and presented here as en example. For information on specific containers for DRAM2 and COMET, please see above.

## 1. Introduction to Singularity Images (Containers)

Singularity containers are lightweight, portable environments that encapsulate an application and its dependencies. They are widely used in scientific computing, bioinformatics, and other fields to ensure reproducibility and simplify software deployment.

## 2. Singularity Containers Used

### COMET Images

- **wetlands-assembly-v7.sif**
  - Description: A wetlands assembly image with a sample name change for coassemblies in 'parse_csv.py'.
  
- **comet-support.sif**
	- Description: Container including samtools, bowtie2, fastqc, multiqc, fastp and sickle.
  
- **wetlands-checkm2-2.sif**
  - Description: Container for running CheckM version 2.

### Gene-Centric Pipeline Images
- **Gene-Centric-Container.sif**
  - Description: Base container for the Gene-Centric pipeline.
  
- **Gene-Centric-Container-v2.sif**
  - Description: Container with 'featureCounts' added.
  
- **Gene-Centric-Container-v3.sif**
  - Description: Container with 'prodigal' included.
  
- **Gene-Centric-Container-v4.sif**
  - Description: Container with additional libraries - 'pandas', 'numpy', and 'openpyxl'.
  
These containers are also available on the server in `/home/opt/singularity_containers/`.

## 3. Usage in Nextflow Script

The Nextflow script uses these Singularity containers to execute different processes within the pipeline. Each process specifies the container to use. These are defined within the `nextflow.config` file.

For example:
```groovy
process {
    withName: COASSEMBLY {
        container = "${params.assembly_container}"
        // Additional process configuration
    }
    // Other processes and their container specifications
}
```

These container variables are defined at the end of the `params{}` section of the `nextflow.config` (**Example**):
```groovy
/* Containers */
	assembly_container = "/home/opt/singularity_containers/wetlands-assembly-v7.sif"
	qc_container = "/home/opt/singularity_containers/wetlands-QC-2.sif"
	qc_bin_container = "/home/opt/singularity_containers/wetlands-Bin-QC-Support-pigz.sif"
	checkm2_container = "/home/opt/singularity_containers/wetlands-checkm2-2.sif"

```

### How to download <New_Wetlands_Pipeline_Name> Singularity Images (Example)
There is a script provided within the <New_Wetlands_Pipeline_Name> GitHub `pull_singularity_containers.sh`:
```bash
singularity pull --arch amd64 library://wrightonlabcsu/default/wetlands-assembly-v7.sif:sha256.7ef413324def70176e6060f0641d91558b88edb26499c4757ee3e3a7ad46806d
mv wetlands-assembly-v7.sif_sha256.7ef413324def70176e6060f0641d91558b88edb26499c4757ee3e3a7ad46806d.sif containers/wetlands-assembly-v7.sif

singularity pull --arch amd64 library://wrightonlabcsu/default/wetlands-checkm2-2.sif:sha256.51b7659d2e6ee7a73f99a0b36c8e0dd4ffa4396c83efd3b14d483b16cb679cb4
mv wetlands-checkm2-2.sif_sha256.51b7659d2e6ee7a73f99a0b36c8e0dd4ffa4396c83efd3b14d483b16cb679cb4.sif containers/wetlands-checkm2-2.sif


```

<New_Wetlands_Pipeline_Name> are told within the installation to run this script:
```bash
sh ./pull_singularity_containers.sh
```

## 4. Changing a Singularity Recipe File

To change a Singularity recipe file, follow these steps:

1. Edit the recipe file with the desired changes.
	1. These recipe files are located on the pipeline's GitHub page within a `containers/` directory.
2. Save the changes.

## 5. Building a Singularity Image

To build a Singularity image from a recipe file:

1. Use the `singularity build` command:
   ```bash
   sudo singularity build <image>.sif <recipe_file>
```

Note:
- You need sudo to run this command
- You cannot use special characters in the names or capital letters when pushing this built image to SyLabs.

## 6. Uploading an Image to Sylabs

Login to SyLabs using the Wrighton Lab's GitHub account:
- Username: `username`
- Password: `password`

To upload a Singularity image to Sylabs, you need an access token. If you don't have one, create it following these steps:

1. Go to [Sylabs Dashboard](https://cloud.sylabs.io/dashboard).
2. Click on "Access Tokens" and create a new access token.

### Access Token (Example, may be expired)

After obtaining an access token, follow these steps to upload the image:

1. Login to SyLabs from the command line:
```bash
	singularity remote login --token <Your_Access_Token>
```

3. Sign the image with your key:
   ```bash
   singularity sign <image>.sif
```

2. Push the image to SyLabs
```bash
singularity push <image>.sif library://wrightonlabcsu/default/<image>.sif
```

## 7. Changing Nextflow Config and pipeline 'pull' script

1. You can change the Nextflow configuration in the `nextflow.config` file. This file contains parameters and settings for the pipeline. Make your desired changes to this file.

For example:
```groovy
/* Containers */
	assembly_container = "/home/opt/singularity_containers/wetlands-assembly-v7.sif"
	qc_container = "/home/opt/singularity_containers/wetlands-QC-2.sif"
	qc_bin_container = "/home/opt/singularity_containers/wetlands-Bin-QC-Support-pigz.sif"
	checkm2_container = "/home/opt/singularity_containers/wetlands-checkm2-2.sif"

```

2. You can change the pipeline's 'pull' script to now pull the updated version of the Singularity image.
	1. First, navigate to the SyLabs library at: https://cloud.sylabs.io/library/wrightonlabcsu and select the image you want to update.
	2. Click on  `Pull Command` to view the command to pull the image.
	   For example:
```bash
singularity pull --arch amd64 library://wrightonlabcsu/default/wetlands-assembly-v7.sif:sha256.f431a320317ff0901ed7bcfb48f212d4d437bcca398bfb0bc9c7e982680d1419
```

3. Using this information, update the rediculously long file name in the 'pull' script.
   For example:
 ```bash
 singularity pull --arch amd64 library://wrightonlabcsu/default/wetlands-assembly-v7.sif:sha256.7ef413324def70176e6060f0641d91558b88edb26499c4757ee3e3a7ad46806d
mv wetlands-assembly-v7.sif_sha256.7ef413324def70176e6060f0641d91558b88edb26499c4757ee3e3a7ad46806d.sif containers/wetlands-assembly-v7.sif

```
