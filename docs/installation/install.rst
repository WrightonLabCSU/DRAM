==============================
Using DRAM2: Quick start guide
==============================

Guide covering the three main commands which encompass a typical DRAM2 run. This includes calling genes from an input MAG(s), annotating the called genes and summarizing these annotations into a distillate.

**Overview:**
   * Brief notes on the DRAM2 command structure
   * Brief notes on DRAM2 output structure
   * Step 1: Calling genes **(dram2 call)**
   * Step 2: Annotating called genes **(dram2 annotate)**
   * Step 3: Summarizing annotations **(dram2 distill)**

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Brief notes on the DRAM2 command structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every DRAM2 use case starts with ``dram2`` followed by a command. All DRAM2 commands and options can be found within the main DRAM2 help menu:

.. code-block:: bash

   dram2 --help

DRAM2 options:

.. code-block:: bash::

  -d, --dram_dir PATH    This is both the location where the output of new
                         DRAM2 actions will go, and also the location where
                         the outputs of past DRAM2 actions and metadata can be
                         found.
  --config_file PATH     Point to a config file that you would like to use for
                         this action specifically. With out of this argument,
                         DRAM2 looks for the config at first in
                         `USER_HOME/.config/dram_config.yaml` then in
                         `/etc/dram2_config.yaml`.
  --version              Show the version and exit.
  -v, --verbose          Verbosity of the logging output, the number of 'v's
                         maps to the default logging level of pythons logging
                         module. The default is 4 aka -vvvv. The mapping is -v
                         = CRITICAL ,-vv = ERROR, -vvv = WARNING, -vvvv =
                         INFO, -vvvvv = DEBUG, -vvvvvv... = NOTSET.
  -t, --threads INTEGER  number of threads to run or the number of processors
                         to use. This command is passed to external tools
                         called by DRAM2, as well as being used in DRAM2
                         itself. Note that increasing the threads may increase
                         the amount of memory required.
  --keep_tmp             Keep all temporary files
  -h, --help             Show this message and exit.

DRAM2 options **must** precede the command you wish to use:

For example, specifying the location of the DRAM2 output and number of threads must precede the ``call`` command and the subsequent ``call`` options:

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> call <call-options>

.. code-block:: bash

   dram2 --help

DRAM2 commands as listed in the help menu:

.. code-block:: bash

   Commands:
     call              Call Genes and Filter FASTAs
     annotate          Annotate Genes with Gene Database
     list_dbs          List available databases
     list_db_sets      List available database sets
     pull_rrna         Pull rRNA Sequences With Barrnap (Not Ready)
     pull_trna         Pull tRNA Sequences With tRNAscan (Not Ready)
     distill           DRAM Distillate
     generate_genbank  Make A DRAM GenBank File (Not Ready)
     merge             Merge DRAM Projects (Gene calls and annotations only)
     strainer          Strain to Genes of Interest (Not Ready)
     neighbors         Pull Genes based on Their Neighborhoods (Not Ready)
     phylotree         Phylogenetic trees to DRAM
     adjectives        Describe Gene Features(Adjectives)
     build_db          Build Your Own Custom DRAM Database
     build_db_list     List the input files you can provied

As you can see from the ``dram2 --help`` output, not all commands listed are functional.

DRAM2 commands, and their corresponding options can be found through their individual ``--help`` menus. 
For example:

.. code-block:: bash

   dram2 call --help

Or

.. code-block:: bash

   dram2 annotate --help

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Brief notes on the DRAM2 output structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As seen above in the help menu output, the DRAM2 option to specify the output directory, ``dram2 -d <path/to/output/directory>``, is used not only to specify the output directory but also specifies the location of previous DRAM2 actions. 

   *It is generally a good idea to keep the same output directory for subsequent DRAM2 actions.*

For example, specifying the same output directory (``-d``) for ``dram2 call`` and ``dram2 annotate``. This is beneficial as this directory will accumulate metadata about your DRAM2 run which expidites subsequent DRAM2 commands using the same input dataset.

^^^^^^^^^^^^^^^^^^^^^
Step 1: Calling genes
^^^^^^^^^^^^^^^^^^^^^

Calling genes results in the creation of a ``genes`` directory populated with a directory for each FASTA input. Each new directory will contain three outputs:

  * ``genes.fna``: FASTA formatted nucleotide sequences of called genes  
  * ``genes.gff``: General Feature Format (GFF3) of called genes
  * ``genes.faa``: FASTA formatted protein sequences of called genes

Bring up the help menu:

.. code-block:: bash

   dram2 call --help

.. code-block:: bash::

   Options:
     -f, --force                        Remove all called genes and information
                                        about them, you will only get the current
                                        set of genes from the command, not the
                                        genes from past runs of call.
     --prodigal_mode [train|meta|single]
                                        Mode of prodigal to use for gene calling.
                                        NOTE: normal or single mode require genomes
                                        which are high quality with low
                                        contamination and long contigs(average
                                        length > 3 Kbp). Read more about this option
                                        in the prodigal wiki:
                                        https://github.com/hyattpd/prodigal/wiki.
     --genes_dir PATH                   The directory to store the genes files to
                                        be used or deleted later. This feature is
                                        beta.
     --prodigal_trans_tables [1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|25]
                                        Prodigal trans tables to use for gene
                                        calling. Read more about this option in
                                        the prodigal wiki:
                                        https://github.com/hyattpd/prodigal/wiki.
     -h, --help                         Show this message and exit.

**Basic usage:**

**Example 1:** Single input FASTA file:

.. code-block:: bash::

   dram2 -d <path/to/output/directory> call <options> /some/path/*.fasta

**Example 2:** For multiple FASTA file inputs in separate directories:

.. code-block:: bash::

   dram2 -d <path/to/output/directory> call <options> /some/path/fasta1.fasta /some/path/fasta2.fasta

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step 2: Annotating called genes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotation of called genes results in the creation of a new directory ``annotated`` which will be populated with an ``raw.tsv`` file.

Bring up the help menu:

.. code-block:: bash

   dram2 annotate --help

.. code-block:: bash::

   Options:
     -s, --use_dbset [metabolism_kegg_set|metabolism_set|adjectives|adjectives_kegg]
     --use_db [camper|cant_hyd|dbcan|fegenie|stats|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref]
                                     Specify exactly which DBs to use. This
                                     argument can be used multiple times, so for
                                     example if you want to annotate with FeGenie
                                     and Camper you would have a command like
                                     `dram2 - o output/dir annotate --use_db
                                     fegenie --use_db camper`, the options
                                     available are in this help.
     --bit_score_threshold INTEGER   The minimum bit score is calculated by a
                                     HMMER or MMseqs search to retain hits.
     --rbh_bit_score_threshold INTEGER
                                     Minimum bit score of reverse best hits to
                                     retain hits.
     --custom_fasta_db_name TEXT     Names of custom databases can be used
                                     multiple times.
     --custom_fasta_db_loc PATH      Location of fastas to annotate against, can
                                     be used multiple times but must match the
                                     number of custom_db_name's.
     --custom_hmm_db_name TEXT       Names of custom hmm databases, can be used
                                     multiple times.
     --custom_hmm_db_loc PATH        Location of HMMs to annotate against, can be
                                     used multiple times but must match number of
                                     custom_hmm_name's
     --custom_hmm_db_cutoffs_loc PATH
                                     Location of file with custom HMM cutoffs and
                                     descriptions, can be used multiple times.
     --tempory_dir PATH              Location of the temporary file where the
                                     annotations will be stored, this file will
                                     still be defeated at the end of the
                                     annotation process if the the tmp flag is
                                     not set.
     -f, --force                     Remove all past annotations and annotate
                                     again.
     -h, --help                      Show this message and exit.

**Basic usage:**

**Example 1:** Annotate using the KEGG database

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> annotate --use_db kegg

**Example 2:** Annotate using multiple databases:

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> annotate --use_db kegg --use_db kegg --use_db kofam --use_db merops

**Example 3:** Annotating with all of the databases which provide entries in the metabolism_summary:

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> annotate --use_dbset metabolism_set

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step 3: Summarizing annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distillation step summarizes the annotated genes within the ``annotated`` directory and generates a new directory ``distill`` which is populated with multiple files:

   * ``genome_stats.tsv``: Genome statistics for all input genomes
   * ``metabolism_summary.xlsx``: Metabolism summary of all input genomes, which gives gene counts of functional and structural genes across a wide variety of metabolisms
   * ``product.tsv``: Coverage of pathways, the coverage of electron transport chain components, and the presence of selected metabolic functions
   * ``product.html``: Interactive heatmap showing coverage of pathways and metabolic functions from the ``product.ts

Bring up the help menu:

.. code-block:: bash

   dram2 distill --help

.. code-block:: bash::

   Options:
     --raw_tsv_path PATH     Location of an raw.tsv. You don't
                                     need to use this option if you are using the
                                     same output_dir for dram with a project config.
                                     If you use this option, you must also use
                                     the force flag to bypass the safeguards that
                                     prevent you from running distill with
                                     insufficient data
     -m, --modules [summarize_metabolism|make_genome_stats|make_product]
                                     What distillate module to run. It can be
                                     time consuming to run all the distillate
                                     module for all projects.
     -f, --force                     Remove skip the normal checks.
     --rrna_path PATH                rRNA output from a dram RNA script. You
                                     don't need to explicitly give this path if
                                     you are using an output_dir from dram with a
                                     project config file. The rRNA run will be
                                     automatically detected if you have a project
                                     config file.
     --trna_path PATH                tRNA output from a dram annotation. You
                                     don't need to explicitly give this path if
                                     you are using an output_dir from dram with a
                                     project meta datafile. The tRNA run will be
                                     automatically detected if you have a project
                                     config file.
     --show_gene_names               Give names of genes instead of counts in
                                     genome metabolism summary. This tool is not
                                     fully supported, and may run into the limits
                                     of Excel. Use with caution.
     --use_db_distilate [camper|cant_hyd|dbcan|fegenie|stats|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref]
                                     Specify exactly which db specific distillate
                                     to use. If you know what you are doing it
                                     may be useful to force the output of the
                                     program. If you already annotated with a
                                     database that has an associated distillate
                                     file eg:methyl and still have your project
                                     config, there should be no need for this
                                     command. If you use this command, you should
                                     have a good idea what you are doing and use
                                     the force command also.
     --custom_summary_form PATH      Custom distillate form to add your own
                                     modules to the metabolism summary. You will
                                     need to read the docs to find the format
                                     that this tsv file must take.
     --genomes_per_product INTEGER   Number of genomes per product.html output.
                                     Decrease value if getting JavaScript Error:
                                     Maximum call stack size exceeded when
                                     viewing product.html in browser. Note that
                                     by default the product html will not be
                                     created if the number of genomes is over
                                     2000. You must pass the make_big_html flag
                                     in order to make that html
     --make_big_html                 It is felt that if the number of genomes is
                                     over 2000 that product may be of limited use
                                     because of the size and the number of html
                                     files that will be made. In order to avoid
                                     the large amount of time it will take to
                                     make these distillates it makes sense to
                                     just make the product html
     -h, --help                      Show this message and exit.

**Basic usage:**

**Example 1:** Basic distillation

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> distill

**Example 2:** Distillation of specific databases.

   *For instance, if you annotated using only KEGG (the same as Example 1 in Annotate):*

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> annotate --use_db kegg


*Then you can specify only the KEGG distillation.*

.. code-block:: bash::

   dram2 -d <path/to/output/directory> -t <#threads> distill --use_db_distilate kegg

