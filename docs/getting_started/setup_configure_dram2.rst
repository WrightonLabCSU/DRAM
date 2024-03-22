Setting Up and Configuring DRAM2
===============================

**This documentaition was imported from DRAM1 it is in an indetermant state of
progress**

**NOTE** If you already have an old release of DRAM2 installed and just want to upgrade you may be able to use your old database.

To install DRAM2 you also must install some dependencies. The easiest way to install both DRAM2 and its dependencies is to use `conda <https://docs.conda.io/en/latest/miniconda.html>`_,  but you can also use manual instructions, or if you are an adventurer you can install a release candidate from this repository .

Conda Installation
------------------

You can downloaded the DRAM2 environment Yaml and install DRAM2 into a new `conda <https://docs.conda.io/en/latest/>`_ environment using this command::

   wget https://raw.githubusercontent.com/WrightonLabCSU/DRAM2/main/environment.yaml
   conda env create -f environment.yaml -n DRAM

If this installation method is used, then all further steps should be run inside the newly created DRAM environment, or with the full path to the executable, use `which` with the environment active to find these, the eg. `which DRAM.py`. This environment can be activated using this command::

   conda activate DRAM

You have now installed DRAM, and are ready to set up the databases.

Pip/Manual Installation

If you do not install via a conda environment, then the dependencies `pandas <https://pandas.pydata.org/>`_, `networkx <https://networkx.github.io/>`_, `scikit-bio <http://scikit-bio.org/>`_, `prodigal <https://github.com/hyattpd/Prodigal>`_, `mmseqs2 <https://github.com/soedinglab/mmseqs2>`_, `hmmer <http://hmmer.org/>`_ and `tRNAscan-SE <http://lowelab.ucsc.edu/tRNAscan-SE/>`_ need to be installed manually. Then you can install DRAM using pip::

   pip install DRAM-bio

You have now installed DRAM, and are ready to set up the databases.

Release Candidate Installation
-------------------------------

The latest version of DRAM is often a release candidate, and these are not pushed to pypi, or Bioconda and so can't be installed with the methods above. You can tell if there is currently a release candidate by reading the `release notes <https://github.com/WrightonLabCSU/DRAM/releases>`_.

To install a potentially unstable release candidate, follow the instructions below. Note the comments within the code sections as there is a context in which commands must be used.

.. code-block:: bash

   # Clone the git repository and move into it
   git clone https://github.com/WrightonLabCSU/DRAM.git
   cd DRAM
   # Install dependencies, this will also install a stable version of DRAM that will then be replaced.
   conda env create --name my_dram_env -f environment.yaml
   conda activate my_dram_env
   # Install pip
   conda install pip3
   pip3 install ./

You have now installed DRAM, and are ready to set up the databases.


Setup Databases
---------------


Download databases
^^^^^^^^^^^

Use the db_builder download command to download databases that are compatible
with your Version of dram.

I Want to Use an Already Setup Databases From DRAM1

If you already installed and set up previous version of dram and want to use your old databases, then you can do it with two steps.

Activate your old DRAM environment, and save your old config::

   conda activate my_old_env
   dram_conf_converter.py my_old_config.txt

I Want to Use an Already Setup Databases From DRAM1
^^^^^^^^^^^

^^^^^^^^^^^

Set up DRAM using the following command::

   prepare_databases --output_dir DRAM_data --kegg_loc kegg.pep

`kegg.pep` is the path to the amino acid FASTA file downloaded from KEGG. This can be any of the gene fasta files that are provided by the KEGG FTP server or a concatenated version of them. `DRAM_data` is the path  to the processed databases used by DRAM. If you already have any of the databases downloaded to your server and don't want to download them again then you can pass them to the `prepare_databases` command by use the `--{db_name}_loc` flags such as `--uniref_loc` and `--viral_loc`.

I don't have access to KEGG
^^^^^^^^^^^

Not a problem. Then use this command::
   DRAM-setup.py prepare_databases --output_dir DRAM_data

Similar to above you can still provide locations of databases you have already downloaded so you don't have to do it
again.

To test that your set up worked use the command `DRAM-setup.py print_config` and the location of all databases provided
will be shown as well as the presence of additional annotation information.

*NOTE:* Setting up DRAM can take a long time (up to 5 hours) and uses a large amount of memory (512 gb) by default. To
use less memory you can use the `--skip_uniref` flag which will reduce memory usage to ~64 gb if you do not provide KEGG
 Genes and 128 gb if you do. Depending on the number of processors which you tell  it to use (using the `--threads`
argument) and the speed of your internet connection. On a less than 5 year old server with 10 processors it takes about
 2 hours to process the data when databases do not need to be downloaded.

Make my database Global
^^^^^^^^^^^

System Requirements
-------------------
DRAM has a large memory burden and is designed to be run on high performance computers. DRAM annotates against a large
variety of databases which must be processed and stored. Setting up DRAM with KEGG Genes and UniRef90 will take up ~500
GB of storage after processing and require ~512 GB of RAM while using KOfam and skipping UniRef90 will mean all
processed databases will take up ~30 GB of disk and will only use ~128 GB of RAM while processing. DRAM annotation
memory usage depends on the databases used. When annotating with UniRef90, around 220 GB of RAM is required. If the KEGG
gene database has been provided and UniRef90 is not used, then memory usage is around 100 GB of RAM. If KOfam is used to
annotate KEGG and UniRef90 is not used, then less than 50 GB of RAM is required. DRAM can be run with any number of
processors on a single node.

Citing DRAM
----------
The DRAM was published in Nucleic Acids Research in 2020 and is available `here <https://academic.oup.com/nar/article/48/16/8883/5884738>`_. If
DRAM helps you out in your research, please cite it.


