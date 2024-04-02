Installation
============

Option 1: Conda Environment
----------------------------

1. Clone the DRAM2 GitHub Repository
2. Download pre-formatted databases
3. `Install Anaconda/Miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
4. `Install Nextflow >= v23.04.2.5870 <https://www.nextflow.io/docs/latest/getstarted.html>`_

Option 2: Singularity Container
-------------------------------

1. Clone the DRAM2 GitHub Repository
2. `Install Singularity >= v3.7.0 <https://docs.sylabs.io/guides/3.0/user-guide/installation.html>`_ (to pull Singularity images from SyLabs).
3. Download Singularity container
4. Download pre-formatted databases
5. `Install Nextflow >= v23.04.2.5870 <https://www.nextflow.io/docs/latest/getstarted.html>`_


General Instructions:
---------------------------

   .. code-block:: bash

      git clone https://github.com/WrightonLabCSU/DRAM2.git
      cd DRAM2
      ./pull_singularity_containers.py
      ./pull_databases_full.py

