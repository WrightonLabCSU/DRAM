Installation
============

1) Clone the DRAM2 GitHub Repository
2) Download Singularity container and pre-formatted databases
3) [Install Nextflow >= v23.04.2.5870](https://www.nextflow.io/docs/latest/getstarted.html)
4) [Install Singularity >= v3.7.0](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) (to pull Singularity images from SyLabs).

#### General Instructions:
--------------------------

.. code-block:: bash

   git clone https://github.com/WrightonLabCSU/DRAM2.git
   cd COMET
   ./pull_singularity_containers.py
   ./pull_databases_full.py

Note for use on HPC:
--------------------

- DRAM2 has SLURM auto submission integrated into the pipeline.
   - To use this feature, ensure you have [SLURM](https://slurm.schedmd.com/quickstart_admin.html) on your cluster.
- **If you do not have SLURM and do not want to use SLURM, use the provided alternative config file: `nextflow-No-SLURM.config`.**
   - **To use this config, you need to add the following to your command: `-c nextflow-No-SLURM.config`.**
