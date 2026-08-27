# DRAM Usage

Here you will find give you basic instructions for running DRAM2
> _The entire list of pipeline configuration parameters can be found here: [Parameters API](params_doc.md). Here we will give some overviews on how to use important parameters and launch DRAM_

## Basic command

Below is an example of common workflow in DRAM2. The user is calling genes from a directory of MAGs, annotating those genes using all available databases, assessing the qualtiy of each MAG. The user is then summarizing their annotations into user-friendly documents, assigning MAG-level metabolic functions, and generating an interactive heatmap of these functions per MAG. This command in being run from the command line in the background. A more detailed description of flags can be found below

`nextflow run BortonWrightonLabs/DRAM --input_fasta [INPUT_FASTA] --outdir [OUTPUT_DIR] --rename --annotate
--anno_dbs all --qc --summarize --traits --visualize -profile singularity -bg` 

## Description of command-line options:

```{note}
All command line options can also be set inside a configuration or params file for set it and forget style parameters instead of needing to be provided each time on the command line. See Nextflow or nf-core documentation for examples.
```

The general command to run DRAM2 is:

`nextflow run BortonWrightonLabs/DRAM [OPTIONS]`

If DRAM2 is installed in a shared location, then the command is:

`nextflow run /path/to/DRAM [OPTIONS]`

If installed in a shared location, replace all example commaands below with `nextflow run /path/to/DRAM [OPTIONS]`

### Launching, Direct or Background

You can launch a Nextflow pipeline, such as DRAM2, directly from the command line. This will output all process information to the current terminal window. If the user is ssh'ed into an HPC and they log out, the run will stop. This can be a fine option if the user is running DRAM2 on a local machine, they launch DRAM2 from a slurm job, or they launch the run in a terminal multiplexer such as `tmux` or `zellij`.

Because the Nextflow scheduler itself takes of very minimal resources, you do not need to launch it from a slurm job, you can launch it in what is called "Background" mode with the `-bg` option. This will allow the user to log out of an ssh session and the DRAM2 run will continue. It will start a process with pid that can be found in a file called `.nextflow.pid` in the launch directory, and if you would like to kill the process you can do so with the command:

`kill "$(cat .nextflow.pid)"`

If you would like to launch DRAM2 in a slurm job, because Nextflow uses very minimal resources, it is suggested to launch it in a small job, such as 1 CPU and 1 GB of RAM. But this job will need to stay alive for the entire duration of the DRAM2 run.

### Important Command-Line Options Explained

`--input_fasta`

This is the location to the input FASTA files. Can be named as such: `*.f*`. File format can be changed with `--fasta_fmt` argument. See [Parameters API](params_doc.md) for more information.

`--outdir`

This is the desired output directory.

`--input_genes`

If the user has already called genes they may use this option to specify the location of a directory containing `*.faa` files. It is key, and is stated in the GitHub documentation, they these files have headers which are unique to a given sample for correct downstream processes. File format can be changed with `--genes_fmt` argument. See [Parameters API](params_doc.md) for more information.


`--annotations`

If the user already has a DRAM2 annotations TSV file, in the correct format, they can provide these using this command-line option.

`--slurm`

This option tells Nextflow to use SLURM as the job scheduler. Additional SLURM options can be specified such as `--partition [PARTITION_NAME]` and `--slurm_node [NODE_NAME]`

### Important Core Nextflow Options

Nextflow provides many command-line options to control how the pipeline is run. Below are some of the most important ones for DRAM2 users. You will also notice that all Nextflow options can be seen by running:

`nextflow run -help`

It is also worth noting that all Nextflow options are specified with a single dash `-`, while all DRAM2-specific options are specified with a double dash `--`.

`-bg`

While the user will still see things being output to the current screen, the user can log out and the DRAM2 run will continue.

`-profile`

This is the Nextflow profile to use. The profile determines how software dependencies are handled and what compute environment settings are used. Common profiles include `singularity`, `docker`, `conda`. The user can also create custom profiles in the `nextflow.config` file.

Additionally, short hand modes exist for common run modes, such as `full_mode`, which will run the entire pipeline (without rename), and with `--anno_dbs all`. See the nextflow.config file on GitHub for the full list of profiles.

`-resume`

If the user has already run this command, or a version of it, Nextflow will look for a `work/` directory, in the current directory, to reuse previous analyses/data. If the user changes command-line options, the pipeline will attempt to resume where these changes were made.

`-with-trace`

This option is a Nextflow-provided option which produces a continuously updated log of DRAM2 processes. This is a good place to check how a run is proceeding and is anything has failed.


#### Pipeline Steps

`--rename`

Rename the headers of the input FASTA files such that they will have a unique prefix based on the FASTA file name. This is optional.

`--annotate`

Annotate called genes. You will need to specify what databases to annotate against with `--anno_dbs [OPTIONS]` or `--use_[DB]`. See the Parameters API documentation for more details, or `--help`.

`--qc`

Run quality control options for DRAM workflow. The QC step collects rRNA and tRNA scans using Barrnap and tRNAscan-SE for the genome_states.tsv output as a baseline. Additional options for QC can be found in the Parameters API documentation or `--help`.

`--summarize`

Distill out annotations from the topic toolkit default (set of predetermined distill topics). You can also specify additional ecosystem toolkits with `--sum_ecos [OPTIONS]` or custom distill sheets. See the Parameters API documentation for more details, or `--help`.

`--visualize`

Create visualizations from the distilled annotations.

`--traits`

Distill out traits from the annotated genes.

### Other Command-Line Options Notes

Other command-line options exist to control specific parameters of the DRAM2 pipeline. These are all described in the [Parameters API](params_doc.md) documentation page, or can be seen by running:

`nextflow run BortonWrightonLabs/DRAM --help`

---

## DRAM2 example commands

Simple run with rename, annotate, QC, summarize, and visualize:

```
nextflow run BortonWrightonLabs/DRAM --input_fasta [INPUT_FASTA] --outdir [OUTPUT_DIR] --rename --annotate
--anno_dbs camper,kegg --qc --summarize --visualize -profile singularity
```

Add resume option and ecosystem summaries:

```
nextflow run BortonWrightonLabs/DRAM --input_fasta [INPUT_FASTA] --outdir [OUTPUT_DIR] --rename --annotate
--anno_dbs camper,kegg --qc --summarize --sum_ecos 'eng_sys,ag' --visualize -profile singularity -resume
```

Run all standard databases and launch on slurm and background:

```
nextflow run BortonWrightonLabs/DRAM --input_fasta [INPUT_FASTA] --outdir [OUTPUT_DIR] --rename --annotate
--anno_dbs all --qc --summarize --sum_ecos 'eng_sys,ag' --visualize -profile singularity -resume --slurm -bg
```

The same as the above command but with full_mode to simplify the command:

```
nextflow run BortonWrightonLabs/DRAM --input_fasta [INPUT_FASTA] --outdir [OUTPUT_DIR] --rename --sum_ecos 'eng_sys,ag' -profile singularity,full_mode -resume --slurm -bg
```

Utilizing a custom nextflow.config file to pass specific parameters (with a custom configuration file, DRAM parameters can be set there and do not need to be specified on the command-line, but Nextflow options still do):

```
nextflow run BortonWrightonLabs/DRAM -c [PATH/TO/NEXTFLOW.CONFIG] -profile singularity -resume -bg
```
