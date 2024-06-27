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

To build a Singularity container, you can use the following command, where `DRAM2-Nextflow-Main-Container-March112024-V4` is the name of your recipe file:

```
sudo singularity build DRAM2-Nextflow-Main-Container-March112024-V4.sif DRAM2-Nextflow-Main-Container-March112024-V4
```

*Note: `sudo` is used here to build the container. While Singularity does have a `--fake-root` option, it is advised to use `sudo`. If a user does not have `sudo` privileges on a machine, it is advised to build them locally, where the user has sudo, and upload them to the machine used for DRAM2 development.*

This command will create a .sif file, which is the Singularity Image Format, containing all the dependencies and software specified in your recipe file. This image can then be used to run processes in your Nextflow pipeline, ensuring a consistent and reproducible environment.

In the `nextflow.config` file on line 296 we define `main_container` as the built `.sif` container:

```
    /* Containers and Environments */
        main_container = "./containers/DRAM2-Nextflow-Main-Container-March252024-V4.sif"
```

Then, `main_container` is subsequently used within the profile definitions:

Example of `singularity_slurm` profile referencing the `main_container`:
```
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
`DRAM2-Nextflow-Main-Container-March252024-V4.sif`

#### DRAM2 visualization container

DRAM2 visualizations (DRAM2 Product) are being implemented as this documentation is being written. Thus, an updated container was created:

`DRAM2-Nextflow-Main-Container-March262024-V5.sif`

To ensure visualizations (DRAM2 Product) can be performed, the `main_container` value must be updated to:

`DRAM2-Nextflow-Main-Container-March262024-V5.sif`

#### Why only one container

This is ideal as images can be large in size. As a software grows, and dependencies become numerous, it can be difficult for all of the dependencies to reside within a single container. However, this can usually be solved by creating Conda environments within a given container. Additionally, installing Conda within a container is typically an easier way to install dependencies within a given container.

Here is an example of installing Conda within a container and performing subsequent installs:

```
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

```
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

```
    /* Containers and Environments */
        main_container = "./containers/DRAM2-Nextflow-Main-Container-March252024-V4.sif"
        main_environment = "./assets/conda/environment.yml"
        quast_environment = "./assets/conda/environment-quast.yml"
```

Here are, the variables `main_environment` and `quast_environment` are defined.

Then, these vairiable are referenced within the definition.

Example Conda profile definition:

```
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

```
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
```
conda = "${params.quast_environment}"
```

While CALL_GENES() uses:
```
conda = "${params.main_environment}"
```

--------

## Future development notes

From development of COMET, it was noticed that relying on Nextflow to build the Singularity containers resulted in conflicts across different HPCs, that is W2 and Riviera. Thus, it was decided, in COMET, to pre-build the containers. Within Nextflow you can, instead of linking the `*.yml` recipe files, provide the location to a pre-build Conda environment. However, this means either the user must build these (As is instructed for COMET) or DRAM2 must provide these.

It is suggested to store the Singularity containers, and pre-built Conda environments if that is chosen, on GLOBUS with the databases. Then a given user can pull everything at once. 