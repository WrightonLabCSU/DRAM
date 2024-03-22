.. DRAM2 documentation master file, created by
   sphinx-quickstart on Wed Nov 30 14:52:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================================
Welcome to the DRAM2 Docs!
=========================================

You are reading the official documentation for DRAM2 (Distilled and Refined Annotation of Metabolism, Version 2) a tool for annotating metagenomic assembled genomes. DRAM2 annotates MAGs  using a set of :ref:`databases<database_info>` including custom user databases. 

DRAM2 is run in three stages. 

   * First an genes are called given a MAG, or set of MAGs, in FASTA format. 
   * Second, an annotation step to assign database identifiers to called genes. 
   * Lastly, a distill step to curate these annotations into useful functional categories. 
   * Additionally, viral contigs are further analyzed during to identify potential AMGs. This is done via assigning an auxiliary score and flags representing the confidence that a gene is both metabolic and viral.

For more detail on how the Science of DRAM2 works please see our DRAM1 `paper <https://academic.oup.com/nar/article/48/16/8883/5884738>`_.

For information on how DRAM is changing, please read the most recent `release notes <https://github.com/WrightonLabCSU/DRAM/releases/latest>`_.

DRAM2 Development Note
----------------------

At this time, DRAM2 is only available internally
The DRAM development team is actively working on DRAM2. We do not anticipate adding any additional functionality to DRAM, i.e. DRAM1. Features requested for DRAM1 will be added to DRAM2, to the best of our ability and as appropriate.




.. toctree::
  :caption: Getting Started
  :name: getting_started
  :hidden:
  :maxdepth: 1

  getting_started/setup_configure_dram2
  getting_started/visualization
  getting_started/basic_use
  getting_started/dram1_to_dram2
  getting_started/update_database
  getting_started/command
  getting_started/example_use_case

.. toctree::
  :caption: Science Guide
  :name: commands
  :hidden:
  :maxdepth: 1

  science/scientific_overview
  science/database

.. toctree::
  :caption: Command Guide
  :name: commands
  :hidden:
  :maxdepth: 1

  commands/dram2
  commands/merger

.. toctree::
  :caption: Developer Documentaition
  :name: dev
  :hidden:
  :maxdepth: 0

  dev/docs
  dev/functions
  dev/making_db_kits
