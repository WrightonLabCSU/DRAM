
# How to run DRAM2 on W2 server

Note: For how to setup DRAM2 on the W2 server refer to the document, How-To-Setup-DRAM2-W2.md.

W2 location to run DRAM2:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF`

-------

## Required files
DRAM2 needs both the `DRAM2.nf` file and the `nextflow.config` files within the directory in which you want to run DRAM2 from. 

DRAM2 also needs the directories which are pulled from the GitHub repository. It is assumed by DRAM2 that the user will run DRAM2 from this location and thus links needed files within the GitHub-pulled directories: `assets` and `modules`.

### Running DRAM2 from a different location

It is suggested that there is a central location on the server where the latest GitHub repository files are located. Then, one must update the locations of the assets, modules and Singularity containers within the files, `DRAM2.nf` and `nextflow.config`.

To dop this, and for more info on setting up DRAM2 in general, see the documentation file, `How-To-Setup-W2-Riviera.md`.

## Description of command-line options:

`nextflow run DRAM2.nf`

This is how to run DRAM2. If you want to immediately run DRAM2 in the background, similar to screen, you can use the `-bg`option:

`nextflow run -bg DRAM2.nf`

While the user will still see things being output to the current screen, the user can log out and the DRAM2 run will continue.

`--input_fasta' 

This is the location to the input FASTA files. Can be named as such: `*.f*`.

`--outdir` 

This is the desired output directory.

`--threads 10`

This sets the threads used for EACH INDIVIDUAL DRAM2 process. This each sample, for each database to be annotated will use 10 CPUs. However, processes such as `HMMSCAN()` can only use 2 CPUs and thus will only use 2 CPUs. This is only a small explanation of a much larger ability to control resource management and more details are in the document, `DRAM2-Computational-Resource-Management.md`. 

`--rename`

Rename the headers of the input FASTA files such that they will have a unique prefix based on the FASTA file name. This is optional.

`--call`

Call genes using Prodigal.

`--annotate`

Annotate called genes.

`--use_camper --use_canthyd`

Annotate the called genes for each input sample using the CAMPER and Cant-Hyd databases.

`--distill_topic default --distill_ecosystem 'eng_sys ag'`

Distill out annotations from the topic toolkit default (set of predetermined distill topics). Also, distill out annotations using the ecosystem toolkits for Engineered Systems and Agriculture. Note, if more than one topic or ecosystem is desired, they much be enclosed in single quotes.

`-profile singularity_slurm`

"singularity" means use the containerized Singularity image (houses all dependencies in exact specified versions) to run DRAM2. 

"slurm" means to use the systems provided SLURM manager to submit, automatically, each DRAM2 process as a SLURM job.

`--no_trees` 

As DRAM2 Trees is still in development, it is suggested to use this command-line option to NOT run DRAM2 Trees.

`-resume`

If the user has already run this command, or a version of it, Nextflow will look for a `work/` directory, in the current directory, to reuse previous analyses/data. If the user changes command-line options, the pipeline will attempt to resume where these changes were made.

`--input_genes`

If the user has already called genes they may use this option to specify the location of a directory containing `*.faa` files. It is key, and is stated in the GitHub documentation, they these files have headers which are unique to a given sample for correct downstream processes.

`--annotations`

If the user already has a DRAM2 annotations TSV file, in the correct format, they can provide these using this command-line option. 

`--queue_size 5`

This option is only present on the last example command. This option allows you to control the total number of SLURM jobs submitted during a given DRAM2 run. This is the easiest way to ensure you limit the consumption of resources on a given server. A queue size of 5 means at a given time, only 5 SLURM jobs will be submitted. Nextflow keeps an internal queue of jobs to submit and the user does not need to be concerned with this. Please see the document, `DRAM2-Computational-Resource-Management.md` for more details on computational resources as DRAM2 can be run at a small scale, queue size = 5, or at an unlimited scale, queue size = 0.

`-with-trace`

This option is a Nextflow-provided option which produces a continuously updated log of DRAM2 processes. This is a good place to check how a run is proceeding and is anything has failed.

`--slurm_node zenith --slurm_queue wrighton-hi`

These are both used to specify the compute node and the SLURM partition. 

-------

## DRAM2 example commands

To be ran in:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF`

Super small - Call, Annotate, Distill:
Has -resume for illustration:
```
nextflow run DRAM2.nf --input_fasta /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/subsampled/super_small/ --outdir DRAM2-super-small-test-2-06252024 --threads 10 --call --rename --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees -with-trace --slurm_node zenith --slurm_queue wrighton-hi -resume 
```


Super small - Call, Annotate, Distill:
NO resume:
```
nextflow run DRAM2.nf --input_fasta /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/subsampled/super_small/ --outdir DRAM2-super-small-test-2-06252024 --threads 10 --call --rename --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees -with-trace --slurm_node zenith --slurm_queue wrighton-hi
```


Super small - Starting from Annotate:
```
nextflow run DRAM2.nf --input_genes DRAM2-super-small-test-2-06252024/Prodigal_v2.6.3/ --outdir DRAM2-super-small-test-2-ANNOTATE-06252024 --threads 10 --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees -with-trace --slurm_node zenith --slurm_queue wrighton-hi
```


Super small - Starting from Distill:
```
nextflow run DRAM2.nf --annotations DRAM2-super-small-test-2-06252024/RAW/raw-annotations.tsv --outdir DRAM2-super-small-test-2-DISTILL-06252024 --threads 10 --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees -with-trace --slurm_node zenith --slurm_queue wrighton-hi
```


Bigger, subsampled - Call, Annotate, Distill:
```
nextflow run DRAM2.nf --input_fasta /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/subsampled/ --outdir DRAM2-super-small-test-2-06252024 --threads 5 --call --rename --annotate --use_kofam --use_kegg --use_dbcan --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 5 -profile singularity_slurm -resume --no_trees -with-trace --slurm_node zenith --slurm_queue wrighton-hi
```


