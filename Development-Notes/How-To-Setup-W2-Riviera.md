# How to Setup DRAM2 on the W2 server and Riviera

*Note 1: This process is laid out within the GitHub page however, when setting up a software for broad use on a server for many users, there are considerations to keep in mind.*

*Note 2: This process is to ideally set up DRAM2 on a server which will serve multiple people. The end goal of this document is to enable a user to setup DRAM2 using this document and then, future users only need to copy the prepared `DRAM2.nf` and `nextflow.config` scripts where they desire to run DRAM2.*

*Note 3: The databases total 687Gb and this process ensures that not every user will download the DRAM2 databases into a given project directory. It is suggested to hold the databases in a central location and modify the `nextflow.config` script to reflect the database location (described below).*

--------

## Setting up DRAM2 on W2

**1) Pull the GitHub repository to a desired location on the server.**


**2) Obtain the databases:**

**Option 1:** (recommended) Replace all occurrences of `./databases/` using the path to the database backup located on W2 to the `nextflow.config` file.

Example for KEGG database at line 140:

Before change:
```
/* Annotation database locations */
        // KEGG
            kegg_db = "./databases/kegg/"
```

After change:
```
/* Annotation database locations */
        // KEGG
            kegg_db = "/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-database-backup-06252024/kegg/"
```

**Option 2:** (not recommended) Copy the databases from the backup location on W2 into the `databases/` directory which was pulled from the GitHub repository. 

Note: This whole database is 687Gb and should not be copied many times!


**3) Set location of `assets/`**

*Note: This only needs to be done in the event other users are going to run DRAM2 from different locations. That is, a given user wants to run DRAM2 from `xyz/` directory and they copy over the `DRAM2.nf` file and `nextflow.config` file which have been set up through this process for them.*

Change all occurenes of `./assets/` in the `nextflow.config` file to have the absolute path to the `assets/` directory (e.g. ${pwd}).

Example for carbon distill topic sheet at line 216:

Before change:
```
        /* Topic */
            distill_carbon_sheet = "./assets/forms/distill_sheets/distill_carbon_Jan252024.tsv"
```

After change:
```
        /* Topic */
            distill_carbon_sheet = "/home/path/to/assets/forms/distill_sheets/distill_carbon_Jan252024.tsv"
```

**4) Set location of `modules/`**

*Note: This only needs to be done in the event other users are going to run DRAM2 from different locations. That is, a given user wants to run DRAM2 from `xyz/` directory and they copy over the `DRAM2.nf` file and `nextflow.config` file which have been set up through this process for them.*

Change all occurrences of `./modules/` within the `DRAM2.nf` file to reflect the location of the `./modules` directory.

**5) Set location of Singularity containers**

Modify line 296 in the `nextflow.config` file to reflect the location of the backup Singularity images on the W2 server.

**Example:**

Before change:
```
    /* Containers and Environments */
        main_container = "./containers/DRAM2-Nextflow-Main-Container-March252024-V4.sif"
```

After change:
```
    /* Containers and Environments */
        main_container = "/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-singularity-container-backups/containers/DRAM2-Nextflow-Main-Container-March252024-V4.sif"
```

**5) Run test help command**

**Run:**
```
nextflow run DRAM2.nf --help
```

```
nextflow run DRAM2.nf --version
```

If there are errors relating to the location of assets, modules or Singularity containers, double check all paths are correct within the `DRAM2.nf` and `nextflow.config` files.

--------

## Setting up DRAM2 on Riviera

Follow the instructions above to pull the GitHub page and update the paths for assets, modules and containers.

Follow the instructions on the GitHub to install Nextflow (can be done with Conda) or, if Riviera has been updated to have a Nextflow module, load the Nextflow module.

The only difference will be that on Riviera, each user will need a copy of the DRAM2 databases and a copy of the Singularity image (Not supported on Riviera, users must use `-profile conda` or `-profile conda_slurm` until Singularity is installed.)

DRAM2 can then be tested the same as above:

**Run:**
```
nextflow run DRAM2.nf --help
```

```
nextflow run DRAM2.nf --version
```

If there are errors relating to the location of assets, modules or Singularity containers, double check all paths are correct within the `DRAM2.nf` and `nextflow.config` files.

*Note: Riviera has specific time limits per job and specific partitions to run on. Thus, a user must use the correct '--time' and `--slurm_partition` parameters in their DRAM2 command-line command.*