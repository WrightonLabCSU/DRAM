===============
DRAM2 Databases
===============

DRAM2 has a wide variety of database for gene annotation. A majority of these databases are described below.

*Some of the databases listed here are in-house databases created by people in the Wrighton Lab and will not be documented here upon release*

A list of DRAM2 databases and database sets can be found within the ``dram2 annotate`` help menu:

.. code-block:: bash::

   dram2 annotate --help

.. code-block:: bash::

   Options:
     -s, --use_dbset [metabolism_kegg_set|metabolism_set|adjectives|adjectives_kegg]
     --use_db [camper|cant_hyd|dbcan|fegenie|stats|kegg|kofam|merops|methyl|heme|pfam|sulfur|uniref]

Alternatively, databases can be listed using the commands: ``dram2 --list_dbs`` or ``dram2 --list_db_sets``.


^^^^^^
CAMPER
^^^^^^

   Usage: ``dram2 annotate --use_db camper``

   A description of CAMPER can be found `here <https://github.com/WrightonLabCSU/CAMPER>`_.

^^^^^^
cant_hyd 
^^^^^^

   Usage: ``dram2 annotate --use_db cant_hyd``

^^^^^^
dbCAN
^^^^^^

   Usage: ``dram2 annotate --use_db dbcan``

   `dbCAN Website <https://bcb.unl.edu/dbCAN/>`_

   Citation: Y. Yin, X. Mao, J. Yang, X. Chen, F. Mao, and Y. Xu, "dbcan: a web resource for automated carbohydrate-active enzyme annotation," Nucleic acids research, vol. 40, no. W1, pp. W445–W451, 2012. `DOI <https://doi.org/10.1093/nar/gks479>`_

^^^^^^
FeGenie
^^^^^^

   Usage: ``dram2 annotate --use_db fegenie``

   `FeGenie GitHub <https://github.com/Arkadiy-Garber/FeGenie>`_

   Citation: Garber AI, Nealson KH, Okamoto A, McAllister SM, Chan CS, Barco RA and Merino N (2020) FeGenie: A Comprehensive Tool for the Identification of Iron Genes and Iron Gene Neighborhoods in Genome and Metagenome Assemblies. Front. Microbiol. 11:37. `DOI <https://doi.org/10.3389/fmicb.2020.00037>`_

^^^^^^^^^^^^
Genome_stats
^^^^^^^^^^^^

   Usage: ``dram2 annotate --use_db stats``

   This kit is used to get genome stats produced by `Prodigal <https://github.com/hyattpd/Prodigal>`_.

^^^^^^
KEGG
^^^^^^

   Usage: ``dram2 annotate --use_db kegg``

   `KEGG Website <https://www.genome.jp/kegg/pathway.html>`_

   Citation:  M. Kanehisa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and M. Tanabe, "Kegg: integrating viruses and cellular organisms," Nucleic acids research, vol. 49, no. D1, pp. D545–D551, 2021. `DOI <https://doi.org/10.1093/nar/gkaa970>`_

^^^^^^
KOfam
^^^^^^

   Usage: ``dram2 annotate --use_db kofam``

   `KOfam Website <https://www.genome.jp/tools/kofamkoala/>`_

    Citation: T. Aramaki, R. Blanc-Mathieu, H. Endo, K. Ohkubo, M. Kanehisa, S. Goto, and H. Ogata, "Kofamkoala: Kegg ortholog assignment based on profile hmm and adaptive score threshold," Bioinformatics, vol. 36, no. 7, pp. 2251–2252, 2020. `DOI <https://doi.org/10.1093/bioinformatics/btz859>`_

^^^^^^
MEROPS
^^^^^^

   Usage: ``dram2 annotate --use_db cant_hyd``

   `MEROPS Website <https://www.ebi.ac.uk/merops/>`_

    Citation: Neil D Rawlings and others, The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a comparison with peptidases in the PANTHER database, Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D624–D632, `DOI <https://doi.org/10.1093/nar/gkx1134>`_

^^^^^^
Methyl
^^^^^^

   Usage: ``dram2 annotate --use_db methyl``

   Methyl is a in-house database mostly made by McKayla Borton.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Heme Regulatory Motifs Counts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   Usage: ``dram2 annotate --use_db heme``

   In-house database.

^^^^^^
Pfam
^^^^^^

   Usage: ``dram2 annotate --use_db pfam``

   `Pfam Website <http://pfam.xfam.org/>`_

    Citation: J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G. A. Salazar, E. L. Sonnhammer, S. C. Tosatto, L. Paladin, S. Raj, L. J. Richardson et al., "Pfam: The protein families database in 2021," Nucleic acids research, vol. 49, no. D1, pp. D412–D419, 2021. `DOI <https://doi.org/10.1093/nar/gkaa913>`_

^^^^^^
Sulfur
^^^^^^

   Usage: ``dram2 annotate --use_db sulfur``

   Generated using the `RefSeq <https://www.ncbi.nlm.nih.gov/refseq/>`_ database.

    Citation: Li W, O'Neill KR, Haft DH, DiCuccio M, Chetvernin V, Badretdin A, Coulouris G, Chitsaz F, Derbyshire MK, Durkin AS, Gonzales NR, Gwadz M, LanczyckiCJ, Song JS, Thanki N, Wang J, Yamashita RA, Yang M, Zheng C, Marchler-Bauer A, Thibaud-Nissen F. RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. Nucleic Acids Res. 2021 Jan 8;49(D1):D1020-D1028. `DOI <https://doi.org/10.1093/nar/gkaa1105>`_

^^^^^^
UniRef
^^^^^^

   Usage: ``dram2 annotate --use_db uniref``

   `UniRef Website <https://www.uniprot.org/help/uniref>`_

    Citation: Y. Wang, Q. Wang, H. Huang, W. Huang, Y. Chen, P. B. McGarvey, C. H. Wu, C. N. Arighi, and U. Consortium, "A crowdsourcing open platform for literature curation in UniProt. PLoS Biol. 2021 Dec 6;19(12):e3001464. `DOI <https://doi.org/10.1371/journal.pbio.3001464>`_

