# Installation Guide

## Requirements

### System Requirements
* Nextflow >= v23.04.2
* One of the following dependency management systems:
  * Conda
  * A Nextflow supported Container Runtime:
    * Apptainer
    * Singularity CE
    * Docker
    * Podman
    * Sarus
* DRAM databases (preformatted and downloaded via Globus)

On many HPC systems, Nextflow, Conda, and some container runtime (usually Singularity or Apptainer) are already installed. If you are using a local system, you will need to install Nextflow and Conda or a container runtime. On an HPC system, these tools are often loaded with modules. Refer to your HPC system's documentation for instructions on how to load modules.

## Installation Steps

### 1. Install Nextflow

If you do not have Nextflow installed, you have two options:

#### Option A: Direct Installation
Follow the instructions from the [Nextflow documentation page](https://www.nextflow.io/docs/stable/install.html).

#### Option B: Conda Installation (Recommended)
Install Nextflow within a conda environment to ensure all dependencies (including Java) are properly installed:
```bash
conda create --name env_nf nextflow -c bioconda
conda activate env_nf
```

### 2. Install Container Runtime (Optional)

If you choose to use a container runtime, make sure you have it installed, working, and available in your PATH. You can verify the installation by running:

```bash
# For Apptainer
apptainer --version

# For Singularity
singularity --version

# For Docker
docker --version

# For Podman
podman --version

# For Sarus
sarus --version
```

If you do not have a container runtime installed, see the [nf-core documentation](https://nf-co.re/docs/usage/getting_started/installation#pipeline-software) for installation instructions.

### 3. Setup DRAM

1. Create a DRAM directory:
```bash
mkdir DRAM
cd DRAM
```

2. Download the DRAM databases:
   - You only need to download the databases you plan to use
   - Access the databases through [Globus](https://app.globus.org/file-manager?origin_id=97ed64b9-dea0-4eb5-a7e0-0b50ac94e889) (UUID: 97ed64b9-dea0-4eb5-a7e0-0b50ac94e889)
   - You will need a Globus account to access and download the data
   - Two options for downloading:
     * Using Globus Web UI: Download databases individually
     * Using [Globus Connect Personal](https://docs.globus.org/globus-connect-personal/install/): Transfer whole folders at once (recommended for bulk transfers)
   - The disk space required depends on which databases you choose to download
   - The download consists of a database folder containing all preformatted databases with the description database

3. Configure DRAM:
   - Download the default configuration file:
   ```bash
   curl -o nextflow.config https://raw.githubusercontent.com/WrightonLabCSU/DRAM/refs/heads/dev/nextflow.config
   ```
   - This file can be customized to change:
     * Database locations
     * Container runtime settings
     * SLURM options
     * Other profile options

### 4. Install and Update DRAM

#### Initial Installation
Install DRAM using Nextflow:
```bash
nextflow pull WrightonLabCSU/DRAM -r dev
```

Note: Once DRAM v2 is out of development, the `-r dev` flag will be removed and the command will be simply:
```bash
nextflow pull WrightonLabCSU/DRAM
```

#### Updating DRAM
To update DRAM, rerun the pull command:
```bash
nextflow pull WrightonLabCSU/DRAM
```

To pull a specific branch or version tag:
```bash
nextflow pull WrightonLabCSU/DRAM -r <branch_or_tag>
```

Note: Only version tags >= 2 are supported. Pre-version 2 releases were before Nextflow was used.

### Dependency Management
DRAM uses Nextflow to handle all dependency management automatically. When you first run DRAM, it will:
1. Download all required dependencies based on your chosen profile (conda, singularity, docker, etc.)
2. Cache these dependencies for future use
3. No additional Python packages or environment setup is needed

## Important Notes

### Installation Location
Nextflow installs all pipeline scripts by default in `$HOME/.nextflow/assets/WrightonLabCSU/DRAM`. If you're running DRAM on a shared system, you may want to install it in a shared directory.

### Custom Installation Location
To install DRAM in a custom location:
1. Clone the DRAM repository to your desired location
2. Run the pipeline by specifying the path:
```bash
nextflow run <path/to/DRAM> <options>
```

### Custom Configuration
To use a custom nextflow.config file:
```bash
nextflow run <path/to/DRAM>/main.nf -c <path/to/custom/nextflow.config> <options>
```

## Verification

To verify your installation, you can run a test command:
```bash
nextflow run DRAM --help
```

This should display the help message with all available options.