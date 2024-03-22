======================
Quick Command Refrence
======================


>> dram2
--------

The root command of DRAM2 protocal which provides infomation to all other commands.

Options:
 * -o, --output_dir PATH  output directory
 * --config_file PATH     Point to a config file that you would like to use.
 * --version              Show the version and exit.
 * -v, --verbose          Verbosity of the logging output, the number of 'v's
 *                        maps to the default logging level of pythons logging
 *                        modual. the default is -vvv. The mapping is 0 =
 *                        CRITICAL ,-v = ERROR, -vv = WARNING, -vvv = INFO,
 *                        -vvvv = DEBUG, -vvvvv... = NOTSET
 * -c, --cores INTEGER    number of processors to use. This will increase the
 *                        amount of memory required
 * --keep_tmp             Keep all temporary files
 * -d, --db_path PATH     Specify a directory Path for any databases that you
 *                        need, but have not setup to live. If you don't
 *                        specify this path and request a databases that has
 *                        not yet been built this will fail. Dram will not use
 *                        files that are in this database if they are not also
 *                        specified in the config file you pass.
 * --help                 Show this message and exit.
 
Commands:
 * call                Call Genes and Filter Fastas...
 * annotate            Annotate Genes with Gene Database ---
 * list_databases      List available database ---
 * pull_rrna           Pull rRNA Sequences With Barrnap (Not Ready) ___
 * pull_trna           Pull tRNA Sequences With tRNAscan (Not Ready) ___
 * distill             DRAM Distillate ---
 * generate_genbank    Make A DRAM GenBank File (Not Ready) ___
 * merger              Merge DRAM Projects (Not Ready) ___
 * strainer            Strain to Genes of Interest (Not Ready) ___
 * neighbors           Pull Genes based on Their Neighborhoods ___ DRAM2...
 * phylotree           Philogenetic trees to DRAM ___
 * adjectives          Describe Gene Features (Adjectives)
 * build_database_set  Build Your Own Custom DRAM Database (Not Ready) ___ 



>> dram2 call
--------

>> dram2 annotate 
--------

>> dram2 list_databases 
--------

>> dram2 pull_rrna
--------

>> dram2 pull_trna
--------

>> dram2 distill
--------

>> dram2 generate_genbank 
--------

>> dram2 merger 
--------

>> dram2 strainer 
--------

>> dram2 neighbors
--------

>> dram2 phylotree
--------

>> dram2 adjectives 
--------

>> dram2 build_database_set 
--------


.. _phylo_trees-how_use:

How do I Use it in DRAM2
===========================


The process of used by phylogentic Tree Explorer in DRAM2 is outlined in the overview figure on the :ref:`phylo_trees-what_how` page.
That figure includes the inputs and outputs of the tool in adition to the process that runs the tool. We will refrance the process coverd in that past

Input Output, and Commands
--------------------------



TODO
---

 - Integrate into dram2 comandline
 - Finalize tests and examples
 - Add all other trees
 - Integrate into new dram2 config (maybe)
 - Alow users to replace the annotations file or make this thing a h5 mofo and kill the user out of that
