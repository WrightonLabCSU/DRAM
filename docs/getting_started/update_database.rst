Updating Databases
===================

Updating the database is supper simple and foremost you will just run::

	dram2 -d <where to put the db> build_db

To put the config in a specific location::

	dram2 -d <where to put the db> build_db -y <yaml_loc>

To update::

	dram2 -d <where to put the db> build_db  -u

To make only a few dbs::

	dram2 -d <where to put the db> build_db  --make_db <db_name>


What about kegg and camper
------------------------

KEGG will always need to be independently downloaded and CAMPER to while its repository is private. The input argument lets you do this.

Download
^^^^^^^

For CAMPER just go to the release page on the github page and download the tag_gz release.

For KEGG use the scripts in this repository, here: https://github.com/WrightonLabCSU/DRAM2/tree/main/extra/kegg_tools
You need to clone the repository or download the zip and copy all the scripts to a
folder on a disk with lots of space. Then you will run them. You can't just get the
script because we have the repo as private.

This is easy, all you need is a username and password. You must run the next
command in screen or tmux, and you will be prompted for your passwords. This will run a
long time and download a bunch of stuff::

	./download_kegg.sh

Now you will have a folder full of files. Don't change directories, don't change the
folder. Run format_kegg.py on the screen or as a big slurm job in the same directory.
You will probably need to activate a dram conda environment, any version of dram, to get
the paths you need. The command is::

	python format_kegg.py

Now you have one pep file "kegg-all-orgs_unique_reheader.pep" ready for DRAM.

Pass File locations
^^^^^^^^^^^^^

To see the arguments and keys for the setting file locations, you can use::
   
	dram2 build_db_list

The command to update camper and KEGG is then::
 
	dram2 -d <output for db> build_db -u --make_db camper --make_db kegg -i <kegg_pep>
	-i camper_tar_gz <path to the tar file>

Do read the help for more info.

