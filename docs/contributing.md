# DRAM: Contribution

Hi there!
Many thanks for taking an interest in improving DRAM.

DRAM is written as a Nextflow Pipeline that wraps around the largely python code that does the DRAM calculations. Nextflow handles a lot of the cross platform deployment, dependency manegement, scheduling, and parallelization. What we need to do to get that out of Nextflow though is write whatever DRAM python code we need that can be interphased with through a CLI so Nextflow can call it directly by passing arguements to it. Then the script should write its outputs to files so Nextflow can pass those arguments to the next step. 

## Nextflow resources 

You can read more about Nextflow [here](https://www.nextflow.io/docs/latest/index.html), 

## Nf-core

DRAM is also following the convention of the nf-core community, which is a a large community bioinformatics Nextflow users that have developed tools and standards for using Nextflow and making it easier to use Nextflow. DRAM is currectly still in the process of transition from the original Nextflow implementation of DRAM to the nf-core style DRAM which will allows more flexibility and easier development.

nf-core docs, tools, and examples can be found [here](https://nf-co.re/docs/)

## Main Structure of DRAM with Nextflow

Pipeline kicks of with `main.nf` which runs some boiler plate initialization steps such as printing the help text if needed, prints your parameters you are using, etc. Then calls the main DRAM pipeline file, which is `workflow/dram.nf`. `workflow/dram.nf` controls the entire flow of the steps of DRAM. Every pipeline step we think of in DRAM runs through this file to some extect. Call Prodigal, Annotate, Distill, Product, etc. When DRAM wants to run Annotate, because Annotate is a large step with many substeps, it then calls `subworkflows/local/annotate.nf`, which then does all of the logic for the Annotate step. Annotate then calls individual processes such as `modules/local/annotate/mmseqs_index.nf` which is then dispatched with the scheduler (SLURM if you are using SLURM, or whatever other scheduling system you are using (local, or whatever)). Those individual processes that get dispatched as individual jobs (or parallel jobs) are what run the DRAM python and other language scripts.

## Pipeline contribution conventions

To make the DRAM code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Adding a new step

1. Write up your main script or script code in whatever language and place it in the bin folder.
    * We need to ensure it is executable for Nextflow, so add a shebang to the top of the file (ex. `#!/usr/bin/env python` for a python file) then run `chmod a+x bin/new_script.py` replacing new_sciprt.py with the actual script.
1. If this is a new process, add a new nextflow script (`.nf`) under `modules/local/whatever_subdir_it_goes` 
1. Define the corresponding input channels and output channels
1. If your scripts in the `bin` folder need external dependencies (such as from conda), then add a `environment.yml` file in the same directory as your process (in the `module/local/..`) in your process script add `conda "${moduleDir}/environment.yml"` above your input channels. (See other process scripts as examples).
1. You will now need to install the wave-cli tool to create containers on demand from conda environments. This is how we specify our containers without having to manage docker files or anything else like that. Wave just tells Nextflow how to pull a container on demand from a conda file. You can install wave-cli [here](https://github.com/seqeralabs/wave-cli). Ideally this will be integrated into our CI/CI system in the future.
1. Once you have wave installed, run `wave --freeze --conda-file modules/local/subdir/environment.yml` which will give you some string that starts something like `community.wave.seqera.io`. Below your `conda "${moduleDir}/environment.yml"`, add `container "YOUR STRING FROM WAVE"`
    * Now users can run DRAM with with `-profile conda` and `-profile singulary/docker/apptainer/etc.` and it will just work without them installing the dependencies
1. Add a computation label to your process, such as `label 'process_low'`, that tells Nextflow how much resources to use. We have defined defaults for what process_low, process_medium, etc. mean, but users can override this in their own configs, allowing more control. All the options can be found in `conf/base.config`
1. Add your process by name to `conf/modules.config` if you want to change where your output files get stored in the user's outdir. 
1. Add a script section to call your DRAM script, passing in whatever CLI argumenets. The outputs from your DRAM script should be the process outputs. You want them to be written directly in the working directory, Nextflow will manage moving them to the users output directy. You can also rename outputs for future Nextflow steps with the `emit:` keyword, and mark some outputs as optional (see other processes).
1. We now need to add this process to the correct part in our pipeline to be called and ran. If it is in an already established step such as Annotate, we might go to `subworkflows/local/annotate.nf` and find the right spot in that code and add it. Though we might need to add a new subworkflow and add that to our `workflows/dram.nf`
1. We also need to add any new parameters to our `nextflow.config` with a default. And then add the equivalent parameter to our `nextflow_schema.json` with help text. The `nextflow_schema.json` is how our `--help` CLI option is populated as well as what our parameters page on our docs is built from. It also allows us to add sanity checks and validation for all relevant parameters.