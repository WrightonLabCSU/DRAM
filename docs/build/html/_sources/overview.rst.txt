Overview
========

DRAM2 (Distilled and Refined Annotation of Metabolism Version 2) is a tool for annotating metagenomic assembled genomes. DRAM2 annotates MAGs using `KEGG <https://www.kegg.jp/>`_ (if provided by the user), `UniRef90 <https://www.uniprot.org/>`_, `PFAM <https://pfam.xfam.org/>`_, `dbCAN <http://bcb.unl.edu/dbCAN2/>`_, `RefSeq viral <https://www.ncbi.nlm.nih.gov/genome/viruses/>`_, `VOGDB <http://vogdb.org/>`_ and the `MEROPS <https://www.ebi.ac.uk/merops/>`_ peptidase database as well as custom user databases. DRAM is run in two stages. First, an annotation step to assign database identifiers to genes, and then a distill step to curate these annotations into useful functional categories. DRAM2 was implemented in `Nextflow <https://www.nextflow.io/>`_ due to its innate scalability on HPCs and containerization, ensuring rigorous reproducibility and version control, thus making it ideally suited for high-performance computing environments. 

For more detail on DRAM and how DRAM works, please see our `paper <https://academic.oup.com/nar/article/48/16/8883/5884738>`_ as well as the `wiki <https://github.com/WrightonLabCSU/DRAM/wiki>`_.

DRAM2 Development Note
----------------------

The DRAM development team is actively working on DRAM2. We do not anticipate adding any additional functionality to DRAM, i.e. DRAM1.

- Future updates will include:

    - Both support for Nextflow + Conda and Nextflow + Singularity (Note: Singularity is not well-supported for MAC.).
    - Pre-formatted annotation and description databases avaiable via `GLOBUS <https://www.globus.org/>`.


Source code is `available on GitHub <https://github.com/WrightonLabCSU/DRAM2/tree/dev>`.


Modules
-------

DRAM2 is organized into 3 main modules. Each module can be run independently of other modules, or strung together in one command.

Call 
^^^^
DRAM2 Call is used to generate called genes from input FASTA files.

    .. code-block:: bash

        nextflow run DRAM2.nf --call 

`Prodigal v2.6.3 <https://github.com/hyattpd/Prodigal>` is utilized here to generate output ``*.fna``, ``*.faa``, and ``*.gff`` files.

Output:

    .. code-block:: bash

        Prodigal_v2.6.3/<*.fna, *.faa, *.gff>

Annotate 
^^^^^^^^

DRAM2 Annotate is used to generate annotations from called genes using both `HMMER v3.3.2 <http://hmmer.org/>` and `MMseqs2 v14.7e284 <https://github.com/soedinglab/MMseqs2>`.

    .. code-block:: bash

        nextflow run DRAM2.nf --annotate <options>

Annotate relies on the user to download pre-formatted databases. Additionally the user can provide pre-formatted custom databases as well.

Output:
    .. code-block:: bash

        RAW/raw-annotations.tsv


Distill 
^^^^^^^

DRAM2 Distill is used to distill down the annotations using expertly-curated distillation input data. Additionally, users may provide custom distill input spreadsheets to distill the annotations.

    .. code-block:: bash

        nextflow run DRAM2.nf .... --distill_topic <option> --distill_ecosystem <option> --distill_custom <option>

Toolkits
--------

DRAM2 toolkits are used to organize commands within the modules based on their scientific cohesion. 


Ingest Toolkit
^^^^^^^^^^^^^^

DRAM2 can ingest a variety of inputs depending on the desired outcome.

Genomic data
~~~~~~~~~~~~

    metaG

    metaT <coming soon>

    metaP <coming soon>

    metaB <coming soon>

Genomic metadata
~~~~~~~~~~~~~~~~

    Taxonomy: 

    Sample, MAG, bin, etc. taxonomy - typically obtained from `GTDB <https://github.com/Ecogenomics/GTDBTk>`

    Quality: 

    MAG or bin quality - typically obtained from `CheckM <https://github.com/Ecogenomics/CheckM>` or `CheckM2 <https://github.com/chklovski/CheckM2>`

    See Input Formats for more details.

Manipulate Toolkit
^^^^^^^^^^^^^^^^^^

Generate GFF and GBK:
~~~~~~~~~~~~~~~~~~~~~

    DRAM2 can generate GFF and GBK output files based on the generated annotations:

        .. code-block:: bash

            nextflow run DRAM2.nf .... --generate_gff --generate_gbk


Strain:
~~~~~~~

    DRAM2 can strain out genes of interest:

    <coming soon>

        .. code-block:: bash

            nextflow run DRAM2.nf .... --strain <options>

Merge and Add annotations:
~~~~~~~~~~~~~~~~~~~~~~~~~~

    DRAM2 can either merge existing DRAM2 annotations together or a user can provide additional annotations to incorporate to a current DRAM2 run.

        .. code-block:: bash

            nextflow run DRAM2.nf --merge_annotations <options>

            nextflow run DRAM2.nf .... --add_annotations <options>


Topic Toolkit
^^^^^^^^^^^^^
    DRAM2 offers expertly-curated distillation of annotations.

        .. code-block:: bash

            nextflow run DRAM2.nf .... --distill_topic <options>

    Current list of available Distill Topics:

    Carbon <description>

    Energy <description>

    Misc. <description>

    Nitrogen <description>

    Transport <description>

    CAMPER <description>

    <Stay tuned for more!>

Ecosystem Toolkit
^^^^^^^^^^^^^^^^^
    DRAM2 offers expertly-curated distillation of annotations.

        .. code-block:: bash

            nextflow run DRAM2.nf .... --distill_ecosystem <options>

    Current list of available Distill Ecosystems:

    Engineered Systems <description>

    Agricultural  <description>

Non-homology Toolkit
^^^^^^^^^^^^^^^^^^^^

Neighborhoods:
~~~~~~~~~~~~~~

    DRAM2 can pull genes based on neighborhoods. <coming soon>

    <description>

Patterns:
~~~~~~~~~~

    DRAM2 can pull genes based on patterns. <coming soon>

    <description>

Phylogenetic Trees:
~~~~~~~~~~~~~~~~~~~

    DRAM2 can construct phylogenetic trees for closely related genes, using precise annotations to explore their evolutionary histories and functional divergences across various organisms. <coming soon>

    <description>

Phenotype Toolkit
^^^^^^^^^^^^^^^^^

Traits:
~~~~~~~

    <description>

Heatmap:
~~~~~~~~

    <description>

Build trait map:
~~~~~~~~~~~~~~~~

    <description>

Multi-Omics Toolkit
^^^^^^^^^^^^^^^^^^^

Paint pathways:
~~~~~~~~~~~~~~~~

    <description>

Analyze pathways:
~~~~~~~~~~~~~~~~~

    <description>

Consensus pathways:
~~~~~~~~~~~~~~~~~~~

    <description>



