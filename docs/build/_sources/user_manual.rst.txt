DRAM2 User Manual
=================

DRAM2 provides various commands for calling genes, annotating called genes, and distilling annotations from input FASTA files. Below are the available commands:

.. code-block:: bash

    nextflow run DRAM2.nf --help (--h)

    nextflow run DRAM2.nf --version (--v)

    Example:

    nextflow run DRAM2.nf --input_fasta ../test_data/ 
        --outdir DRAM2-test-data-Feb012024/ 
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

Bring up help menu(s)
----------------------------
.. code-block:: bash

    nextflow run DRAM2.nf --help (--h)


Bring up versions menu
----------------------------
.. code-block:: bash


    nextflow run DRAM2.nf --version (--v)


----------------------------------------------------------------------------------------------------------------------------

Input data
----------------------------

Starting from scratch: Provide DRAM2 with input FASTA files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FASTA files:
^^^^^^^^^^^^

- FASTA files must be located in a single directory
- A given FASTA file represents one sample/MAG/isolate
- FASTA files have one of the following extensions: .fa, .fasta, .fna
- FASTA headers are either renamed with respect to the individual sample or you will plan to use the `--rename` option

(Optional) Additional files:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sample taxonomy:
^^^^^^^^^^^^^^^^

- Sample taxonomy input (``--taxa <path/to/file>``) is modeled after `GTDB-Tk output <https://ecogenomics.github.io/GTDBTk/files/summary.tsv.html#files-summary-tsv>`_
- Minimum requirements: 2 columns named, "user_genome" (sample names matching input file names) and "classification" (taxonomy)
- Taxonomy information will be added to the output `raw-annotations.tsv` and `distillate.xlsx`

+--------------+----------------+
| user_genome  | classification |
+==============+================+
| sample1      | taxonomy1      |
+--------------+----------------+
| sample2      | taxonomy2      |
+--------------+----------------+
| sample3      | taxonomy3      |
+--------------+----------------+


**Note:** DRAM2 visualization may rely on GTDB-Tk-formatted taxonomy and it is in your best interest to follow this format. 

Sample Bin/MAG quality:
^^^^^^^^^^^^^^^^^^^^^^^

- Sample MAG/Bin quality input (``--bin_quality <path/to/file>``) is modeled after `CheckM/CheckM2 output <https://github.com/chklovski/CheckM2>`_
- Minimum requirements: 3 columns named, "Bin Id" (sample names matching input file names), "Completeness" and "Contamination"
- Bin quality information will be added to the output `raw-annotations.tsv` and `distillate.xlsx`

+---------+--------------+--------------+
| Bin Id  | Completeness | Contamination|
+=========+==============+==============+
| sample1 | 95%          | 2%           |
+---------+--------------+--------------+
| sample2 | 88%          | 5%           |
+---------+--------------+--------------+
| sample3 | 92%          | 1.5%         |
+---------+--------------+--------------+

Starting with Annotation: Provide DRAM2 with input called genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Input called genes (proteins), ideally from the latest version of Prodigal
- Input called genes must have the extension: .faa
- Input called gene headers contain the sample-name prefix (this is crucial to ensure annotations are sample-specific) 


Example:
Sample01.faa has the following First entry:

    .. code-block:: bash

        >Sample01_k95_149252_1 # 1 # 513 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.392
        VVNQLANNRQGTQSAKTRSEVSGGGRKPWRQKGTGHARQGSTRSPQWTGGGVVFAPKPRD

**Note:** DRAM2 uses rename.sh from `BBTools <https://jgi.doe.gov/data-and-tools/bbtools/>`_

Starting with Distill: Provide DRAM2 with previous annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Input annotations (TSV) file must be generated via DRAM2
- (Optional) rRNAs (``--rrnas <path/to/file>``) identified with `Barrnap <https://github.com/tseemann/barrnap>`_
- (Optional) tRNAs  (``--trnas <path/to/file>``) identified with `tRNAscan-SE  <http://lowelab.ucsc.edu/tRNAscan-SE/>`_
- tRNA and rRNA information is included in the output `distillate.xlsx`

----------------------------------------------------------------------------------------------------------------------------

General Command-line Options
----------------------------

The usage of DRAM2 involves various commands and options, as outlined below:

- To perform a combined operation of renaming, calling, annotating, and distilling, use:

    .. code-block:: bash

        nextflow run DRAM2.nf --rename --call --annotate --use_<database(s)> --distill_topic <distillate(s)>

- To call genes using input FASTA files and optionally rename the FASTA headers, use:

    .. code-block:: bash

        nextflow run DRAM2.nf --call --rename --input_fasta_dir <path/to/fasta/directory/>

- To annotate called genes using input FASTA files, use:

    .. code-block:: bash

        nextflow run DRAM2.nf --annotate --input_genes <path/to/called/genes/directory>

- To distill annotations using input annotation files, use:

    .. code-block:: bash

        nextflow run DRAM2.nf --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>

- For a more realistic example with multiple options, see:

    .. code-block:: bash

        nextflow run DRAM2.nf 
            --input_fasta ../test_data/ 
            --outdir DRAM2-test-data-Feb012024/ 
            --call --rename 
            --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur 
            --add_annotations ../test-data/old-DRAM2-annotations.tsv
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
            -profile conda_slurm
            -slurm_node alpha

----------------------------------------------------------------------------------------------------------------------------

Call Command-line Options
----------------------------

Call description
~~~~~~~~~~~~~~~~~
The purpose of DRAM2 --call is to call genes on input FASTA files.

Usage
~~~~~
To call genes using input FASTA files, use the following command:
.. code-block:: bash

    nextflow run DRAM2.nf --call --input_fasta_dir <path/to/fasta/directory/> --outdir <path/to/output/directory/> --threads <threads>

REQUIRED DRAM2 profile options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``-profile``: STRING <conda, conda_slurm, singularity, singularity_conda>

    Runs DRAM2 either using Conda (must be installed) or Singularity (must be installed).
    Runs DRAM2 with no scheduling or scheduling via SLURM.
    See SLURM options in the full help menu.

Call options
~~~~~~~~~~~~
- ``--rename``: Rename FASTA headers based on file name.

    Why? DRAM2 output is focused on scaffolds/contigs with respect to each provided input sample.
    Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguishable.
    *If you have already renamed your FASTA headers, do not include '--call'*


- ``--prodigal_mode``: STRING <single|meta>

    Default: 'single'

- ``--prodigal_tras_table``: <1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25>

    Specify a translation table to use (default: '1').

- ``--min_contig_len``: NUMBER <number in base pairs>

    Default: '2500'

----------------------------------------------------------------------------------------------------------------------------

Annotate Command-line Options
-----------------------------

Annotate description
~~~~~~~~~~~~~~~~~~~~
The purpose of DRAM2 '--annotate' is to annotate called genes on input (nucleotide) FASTA (fa*) files.

Usage
~~~~~
To annotate called genes, use one of the following commands:
.. code-block:: bash

    Annotate called genes using input called genes and the KOFAM database:
    nextflow run DRAM2.nf --annotate --input_genes <path/to/called/genes/directory> --use_kofam

    Annotate called genes using input fasta files and the KOFAM database:
    nextflow run DRAM2.nf --annotate --input_fasta <path/to/called/genes/directory> --use_kofam

REQUIRED DRAM2 profile options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``-profile``: `<conda, conda_slurm, singularity, singularity_conda>`

    Runs DRAM2 either using Conda (must be installed) or Singularity (must be installed).
    Runs DRAM2 with no scheduling or scheduling via SLURM.
    *See SLURM options in the full help menu.*



Annotate options
~~~~~~~~~~~~~~~~

- ``--use_<db-name>``: STRING <camper|cant_hyd|dbcan|fegenie|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref>

    Specify databases to use. Can use more than one. Can be used in combination with --use_dbset.

- ``--use_dbset``: STRING <metabolism_kegg_set|metabolism_set|adjectives_kegg_set|adjectives_set>

    metabolism_kegg_set = kegg, dbcan, merops, pfam, heme
    metabolism_set        = kofam, dbcan, merops, pfam, heme
    adjectives_kegg_set = kegg, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
    adjectives_set        = kofam, dbcan, merops, pfam, heme, sulfur, camper, methyl, fegenie
    *Only one set can be used. Can be used in combination with --use_[db-name]*

- ``--add_annotations``: PATH <path/to/old-annoations.tsv>

    Used to add in old annotations to the current run. (See example for format.)

- ``--generate_gff``: OPTION

    Will generate an output GFF for each sample based on the raw-annotations.tsv.

- ``--generate_gbk``: OPTION

    Will generate an output GBK for each sample based on the raw-annotations.tsv.

----------------------------------------------------------------------------------------------------------------------------

Distill Command-line Options
----------------------------

Distill description
~~~~~~~~~~~~~~~~~~~~
The purpose of DRAM2 --distill is to distill down annotations based on curated distillation summary form(s). 
User's may also provide a custom distillate via --distill_custom <path/to/file> (TSV forms).
Distill can be run independently of --call and --annotate; however, annotations must be provided (--annotations <path/to/annotations.tsv>). 
Optional tRNA, rRNA, and bin quality may also be provided.

Usage
~~~~~
To distill annotations, use the following command:
.. code-block:: bash

    nextflow run DRAM2.nf --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv> --outdir <path/to/output/directory/> --threads <threads>

*Important: if more than one topic or ecosystem is included, they must be enclosed in single quotes. Example: --distill_topic 'carbon transport'*

Example
~~~~~~~
Call and Annotate genes using input fastas and KOFAM database. Distill using carbon topic and AG ecosystem:
.. code-block:: bash

    nextflow run DRAM2.nf --input_fasta_dir <path/to/fasta/directory/> --outdir <path/to/output/directory/> --call --annotate --distill_topic carbon --distill_ecosystem ag --threads <threads> --use_kofam

REQUIRED DRAM2 profile options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``-profile``: STRING <conda, conda_slurm, singularity, singularity_conda>

    Runs DRAM2 either using Conda (must be installed) or Singularity (must be installed).
    Runs DRAM2 with no scheduling or scheduling via `SLURM <https://slurm.schedmd.com/documentation.html>`_.
    See SLURM options in the full help menu.


Distill options
~~~~~~~~~~~~~~~~
- ``--annotations``: PATH <path/to/annotations.tsv>

    Required if you are running distill without --call and --annotate.

- ``--rrnas``: PATH <path/to/rRNA.tsv> (See example for format.)

    rRNA information will be included in distill output.

- ``--trnas``: PATH <path/to/tRNA.tsv> (See example for format.)

    tRNA information will be included in distill output.

- ``--bin_quality``: PATH <path/to/bin-quality.tsv> (See example for format.)

    CheckM and CheckM2 compatible.

- ``--taxa``: PATH <path/to/bin-taxonomy.tsv>

    Compatible with GTDB. (See example for format.)

- ``--distill_topic``: STRING <carbon|energy|misc|nitrogen|transport> OR <default = carbon, energy, misc, nitrogen, transport>

    If more than one topic included, they must be enclosed in single quotes.

- ``--distill_ecosystem``: STRING <eng_sys|ag>

    If more than one ecosystem included, they must be enclosed in single quotes.

- ``--distill_custom``: STRING <path/to/custom_distillate.tsv> (See example for format and options.)

    As of now, only one custom distillate may be included.


----------------------------------------------------------------------------------------------------------------------------

Select a Profile which suits your compute setup
------------------------------------------------

DRAM2 utilizes either Conda or Singularity for dependency management and the user MUST choose one of the following options on execution of any DRAM2 command

*The Nextflow profile option is used (`-profile`) - yes! a single hyphen!*

1. 

    .. code-block:: bash
    
        -profile conda


    This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

2. 

    .. code-block:: bash
    
        -profile conda_slurm

    This option will submit each individual DRAM2 process as its own SLURM job. (See Wiki Resource Management for details).
    This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

3. 

    .. code-block:: bash
    
        -profile singularity

    This option relies on the local systems Singularity. Nextflow will create its own Conda environments to run in. 

4. 

    .. code-block:: bash
    
        -profile singularity_slurm

    This option will submit each individual DRAM2 process as its own SLURM job.
    This option relies on the local systems Singularity to run the downloaded Singularity container.  

Which is Better?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Conda Environments
^^^^^^^^^^^^^^^^^^

**Pros:**

- Beginner-Friendly: Easy to install and use, making it accessible for newcomers.
- Reproducibility: Efficient management of environments facilitates reproducibility.

**Cons:**

- Generally slower than using Singularity containers (will have metrics in the future).
- Dependency Conflicts: Dependency resolution can be slow and may lead to conflicts.
- Limited Portability: System dependencies may introduce variability, affecting portability.
- System Variability: Reliance on the host system's architecture and libraries can cause variability between systems.

Singularity Containers
^^^^^^^^^^^^^^^^^^^^^^

**Pros:**

- Generally faster than using Conda environments (will have metrics in the future).
- Consistent Environments: Ensures consistent runtime environments, enhancing reproducibility.
- HPC Ideal: Perfect for high-performance computing (HPC) environments without the need for root access.
- Isolation: Offers isolation from the host system, minimizing conflicts.
- Wide Portability: Containers are portable across any Linux system with Singularity.

**Cons:**

- Installation Complexity: Can be trickier to install compared to Conda environments.
- Storage Space: May consume more storage space.

**Summary**

Conda is recommended for its ease of use and versatility across different programming languages.
Singularity excels in ensuring reproducibility and compatibility in high-performance computing environments.


----------------------------------------------------------------------------------------------------------------------------

DRAM2 Databases
----------------------------

Unlike DRAM1 databases, DRAM2 databases will be pre-formatted and hosted online. Users of DRAM2 will need to:

1. Decide which databases suit their needs.
2. Download DRAM2 databases via the provided ``pull_databases_*.py`` scripts.


**These databases can be quite large, so it's important to review the options below.**

*These databases rely on an SQL database of database descriptions, provided in three different sizes based on the user's needs.*

All databases DRAM2 accommodates:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`KEGG <https://www.genome.jp/kegg/>`_

- Kyoto Encyclopedia of Genes and Genomes.
- (140G)

`dbCAN <http://bcb.unl.edu/dbCAN2/>`_

- A database for automated carbohydrate-active enzyme annotation.
- (202M)

`Kofam <https://www.genome.jp/tools/kofamkoala/>`_

- Customized HMM database of KEGG Orthologs (KOs).
- (14G)

`MEROPS <https://www.ebi.ac.uk/merops/>`_

- A database of proteolytic enzymes and their substrates.
- (3.6G)

`Viral <https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239>`_
- RefSeq viral database.
- (1.6G)

`CAMPER <https://github.com/WrightonLabCSU/CAMPER>`_

- Curated Annotations for Microbial (Poly)phenol Enzymes and Reactions.
- (846M)

`CANT-HYD <https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation>`_

- Curated database of phylogeny-derived hidden markov models for annotation of marker genes involved in hydrocarbon degradation.
- (877M)

`FeGenie <https://github.com/Arkadiy-Garber/FeGenie>`_

- HMM-based identification and categorization of iron genes and iron gene operons in genomes and metagenome assemblies.
- (6.6M)

`Sulfur <url_to_sulfur_placeholder>`_

- Custom Sulfur database.
- (1.7M)

`Methyl <url_to_methyl_placeholder>`_

- Hidden Markov models (HMMs) based on genes related to iron acquisition, storage, and reduction/oxidation in Bacteria and Archaea.
- (52K)

`UniRef <https://www.uniprot.org/help/uniref>`_

- A comprehensive and non-redundant database of protein sequences.
- (477G)

`Pfam <https://pfam.xfam.org/>`_

- A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).
- (8.8G)

`VOGDB <url_to_vogdb_placeholder>`_

- Placeholder description.
- (4.5G)


Downloadable DRAM2 Database Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Big Set**

- Includes all medium databases + UniRef
- Excludes KEGG*

    .. code-block:: bash

        ./pull_databases_full.py

OR

Follow these instructions to pull manually via `GLOBUS <https://www.globus.org/>`_.

**Routine Set**

- Includes: dbCAN, Kofam, MEROPS, Viral, CAMPER, CANT-HYD, FeGenie, Sulfur, Methyl, Pfam, VOGDB
- Excludes KEGG*
- Excludes UniRef

    .. code-block:: bash

        ./pull_databases_routine.py

OR

Follow these instructions to pull manually via `GLOBUS <https://www.globus.org/>`_.

**Minimal Set**

- Includes: <TBD>

    .. code-block:: bash

        ./pull_databases_minimal.py

OR

Follow these instructions to pull manually via `GLOBUS <https://www.globus.org/>`_.

----------------------------------------------------------------------------------------------------------------------------

SLURM Options
-------------

DRAM2 can utilize `SLURM <https://slurm.schedmd.com/documentation.html>`_ for scheduling jobs.

SLURM is used when either of these ``--profile`` options are used: `conda_slurm` or `singularity_slurm`

- Each DRAM2 process is submitted as an individual SLURM job
- This means DRAM2 can scale horizontally on an HPC

SLURM command-line options
~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify SLURM job duration:

    ``--time <time in SLURM format>`` 

    Example: ``--time 10h``

    Example: ``--time 7d``

Specify a specific node to compute on:

    ``--slurm_node <node_name>``

    Example: ``--slurm_node alpha``

Specify a specific partition/queue to compute within:

    ``--slurm_queue <node_name>`` 

    Example: ``--slurm_queue smith``

    Example: ``--slurm_queue 'smith-hi,smith-low'``


**Note:** SLURM can be tricky because administrators do not set up SLURM the same across machines.
You may need to modify a given profile `.config` file to get DRAM2 with SLURM to work on your HPC.

**For Example**, a HPC our group computes on does not allow for memory specificaton. Thus, to use DRAM2 with SLURM on this HPC, we comment out ("//") the ``memory=`` lines in the config files.

The config files are located in:
``/assets/conda/conda_slurm.config`` and ``/assets/singularity/singularity_slurm.config``

----------------------------------------------------------------------------------------------------------------------------

Software Used
-------------

- **BBTools** `v39.01 <https://jgi.doe.gov/data-and-tools/bbtools/>`_
- **Bowtie2** `v2.5.1 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_
- **Prodigal** `v2.6.3 <https://github.com/hyattpd/Prodigal>`_
- **Python** `v3.10 <https://www.python.org/downloads/release/python-3100/>`_
- **Pandas** `v1.5.2 <https://pandas.pydata.org/pandas-docs/version/1.5.2/>`_
- **Pytest** `v7.2.0 <https://docs.pytest.org/en/7.2.x/>`_
- **Scikit-bio** `v0.5.7 <http://scikit-bio.org/>`_
- **MMseqs2** `v14.7e284 <https://github.com/soedinglab/MMseqs2>`_
- **HMMER** `v3.3.2 <http://hmmer.org/>`_
- **SciPy** `v1.8.1 <https://www.scipy.org/>`_
- **SQLAlchemy** `v1.4.46 <https://www.sqlalchemy.org/>`_
- **Barrnap** `v0.9 <https://github.com/tseemann/barrnap>`_
- **Altair** `v4.2.0 <https://altair-viz.github.io/>`_
- **OpenPyXL** `v3.0.10 <https://openpyxl.readthedocs.io/en/stable/>`_
- **NetworkX** `v2.8.8 <https://networkx.org/>`_
- **Ruby** `v3.1.2 <https://www.ruby-lang.org/en/downloads/>`_
- **GNU Parallel** `v20221122 <https://www.gnu.org/software/parallel/>`_
- **tRNAscan-SE** `v2.0.12 <http://lowelab.ucsc.edu/tRNAscan-SE/>`_
- **Samtools** `v1.17 <http://www.htslib.org/>`_
- **CD-HIT** `v4.6 <http://weizhong-lab.ucsd.edu/cd-hit/>`_
- **CoverM** `v0.6.1 <https://github.com/wwood/CoverM>`_
- **Subread** `v2.0.6 <http://subread.sourceforge.net/>`_
- **XlsxWriter** `v3.1.6 <https://xlsxwriter.readthedocs.io/>`_
- **Numpy** `v1.26.0 <https://numpy.org/>`_
