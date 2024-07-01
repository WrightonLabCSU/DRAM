# DRAM2 Computational Resource Management

Implementation of DRAM2 in Nextflow enables a vast array of computational resource management options.

This document will breakdown 1) What Nextflow offers for computational control and what is meant by scaling horizontally, 2) What has been implemented in DRAM2 in terms of managing the scale of a DRAM2 run, 3) Some examples with explanations of various command-line options, and 4) Suggestions on improvements to the current implementation.

------


## What does Nextflow offer

Nextflow offers **many** options for control and reporting of computational resources. Here we are mainly talking about CPUs and memory (RAM) however, this could also include time and GPU usage.

### Nextflow and horizontal scaling

This concept is easiest to grasp in terms of DRAM1 or any other linear pipeline. 

When a user runs a linear pipeline each sample in processed 1-by-1. Thus 5 samples would have their genes called via Prodigal, one at a time, and then these called genes are annotated one at a time. Thus, we can think of 5 samples taking 5x the time one sample would take. While this is by no means wrong, we are losing out on efficiency. For example, when the user states in DRAM1 that they would like to use 20 CPUs, each sample gene calling or annotation is allocated 20 CPUs. Most likely the user ran this in a bash script and submitted through SLURM with 20 CPUs. Where the inefficiency lies is some dependencies can not use 20 CPUs, for example: hmmscan can only use 2 CPUs. Thus, while hmmscan is running, there are 18 CPUs which are reserved but not being used.

Now lets think of this example taking into consideration Nextflows ability to scale horizontally. 

If the same user provides the same 5 samples, DRAM2 with Nextflow will submit a separate job (a SLURM job if the user chose to use SLURM) for EACH samples gene calling and called gene annotation. Thus, 5 gene calling jobs will be submitted at the same time and 5 annotation jobs will be submitted at the same time. Therefore, depending on the CPUs and memory avaiable on the HPC/Server, all of these jobs could occur at the same time. Keep in ming this is an example and gene calling and gene annotation cannot occur at the same time because one is dependent on the other. But for example 1/5 gene calls could be done and while those are running, 1/5 gene annotations would start running in parallel. 

How is this more efficient? hmmscan is a great example: instead of a user specifying they want 20 CPUs for a WHOLE DRAM2 run, the user specifies the max CPUs a single job can use, `--threads 20`, but DRAM2 will only allocate 2 CPUs to hmmscan, ensuring CPUs are not reserved and not used.

Extending this to the extreme with an example:

A user wants to annotate 100 FASTA files with 10 unique databases (e.g. KEGG, CAMPER, etc.) with 10 CPUs. This means there is the potential to be annotating 100 Fasta files * 10 databases * 10 CPUs = 1000 jobs occurring at the same time = 10,000 CPUs being used.

While this is an extreme example, it is used purely to illustrate how DRAM2, implemented in Nextflow, can vastly scale and manage resources.

In the following sections there are implemented command-line arguments and options to control this scalability so DRAM2 can be run on a laptop at a small scale, on a server like W2 with moderate resources, or Alpine with potentially vast resources.

### Nextflow Report Generation

Nextflow offers an array of reports which can aid in optimizing a DRAM2 run or the software itself. 

`-with-report`

This report has computational resource statistics and metrics for each individual DRAM2 process for a given DRAM2 run. This report is generated at the end of a DRAM2 run, successful or not.

`-with-timeline`

This report gives a timeline-based view, with statistics and metrics, for a given DRAM2 run. This report is generated at the end of a DRAM2 run, successful or not.

`-with-trace`

This is a minimal report which is continually updated during a DRAM2 run. This report details the working directory locations for each individual DRAM2 process as well as some resource usage statistics and metrics.
**It is suggested to use this option on EVERY DRAM2 run.**


### Nextflow CPU and memory allocation

Nextflow submits each process as a separate job to the computer. This is be done with or without SLURM.

Each process is defined within the DRAM2 `nextflow.config` file (or the location of the config is linked to within the `nextflow.config` file).

Here is an example of the HMM_SEARCH() process which is used to query called genes against a given HMM database:

(Taken from `./assets/singularity/singularity_slurm.config`)
```bash
    withName: HMM_SEARCH {
        publishDir = [
            path: "${params.outdir}/HMM_search",
            mode: 'copy',
        ]
        
        executor="slurm"
        cpus = '2'
        time = "${params.time}"
        memory = "${params.sml_job}"
        queue = "${params.slurm_queue}"
        clusterOptions = params.slurm_node ? "--nodelist=${params.slurm_node} --job-name=DRAM2-hmm-search" : "--job-name=DRAM2-hmm-search"

        maxForks = "${params.max_forks_hmm}"
    }
```

Let's break it down:

**`publishDir`**

This is the location where files from this process will be output relative to the user's specified output directory `--outdir`.

**`executor`**

This is where we state we want to use SLURM and submit each of these processes as an individual SLURM job. For example, if we look at the same process definition in `./assets/conda/conda_slurm.config`, we will see `executor="conda"`.

**`cpus`**

This is the number of CPUs this process can use. For example, in processes which can reliably use any number of CPUs, like MMSeqs2 query searched, we would have the user specified number of CPUs (`--threads 20`) `cpus = ${params.threads}`.

**`time`**

This is the amount of time a given hmm_scan process can take. On HPCs/servers which do not have time limits, a default of `148hr` is used. On HPCs which have strict time limits, for example the CSU Rivera server, the user needs to account for this on the command line: `--time 12h`.

**`memory`**

Similar to CPU usage. DRAM2 has defaults set for each process based on prior knowledge of each dependent software. 

Within the `nextflow.config` there are defaults set for each process based on prior knowledge of each dependent software:

```bash
    /* Memory allocation */
        xsml_job = '200MB'
        sml_job = '12GB'
        med_job = '50GB'
        lrg_job = '200GB'
        xlrg_job = '400GB'
```

**`clusterOptions`**

Cluster options is essentially any supported SLURM parameter which can be set. Here there is `--nodelist` and `--job-name`. Additionally there is `slurm_queue` to set a particular SLURM partition. The `--job-name` is the name which will appear in a SLURM queue if a user checks `squeue`. 

**`maxForks`**

This is where Nextflows ability to scale horizontally comes into to play. Back to the example above where the user provided 5 input FASTA file. If `maxForks` for the HMM_SCAN() process is set to 5, then all 5 samples can occur at the same time. However, if `maxForks` is set to 1, then only 1 of the samples can be run through HMM_SCAN() at a time. 

Within the `nextflow.config` there are defaults set for each process based on prior knowledge of each dependent software and its computational demands:

```bash
    /* Max forks Options */
        max_forks_single_cpu  = 10
        max_forks_user_cpu = 2
        max_forks_hmm = 5
```

Essentially, if a process can only use 1 CPU, then we will at most run 10 of these processes at a time. If a process can use a varibale number of CPUs, based on the user's command-line options (`--threads 10`) then only 2 of these will occur at a time. Here we have the HMM_SCAN() process hard-coded to run 5 of these at a given time. 

Thus, it is crucial a user considers how many samples they have and what computational resources they have available. 

*Note: this may seem like too much to think about just to annotate some FASTA file. It is! That is why there is one additional command-line option provided through Nextflow to control the total number of processes which can occur at a given time (--queue_size ) and is described below.

------

## Implentation of controllable resources in DRAM2

### Threads/CPU options

For a given DRAM2 run the user has the basic ability to set the max number of CPUs a given job will be able to use: `--threads 10`.

**Note: this does not mean, "For a total DRAM2 run I as a user will only use 10 CPUs at a time". This correctly means, "For each process in DRAM2 (calling genes or annotating called genes) will use at most 10 CPUs.

### Memory options

Within the `nextflow.config` there are defaults set for each process based on prior knowledge of each dependent software:

```bash
    /* Memory allocation */
        xsml_job = '200MB'
        sml_job = '12GB'
        med_job = '50GB'
        lrg_job = '200GB'
        xlrg_job = '400GB'
```

These can be modified within the `nextflow.config` file itself to limit resources, or can be added at the command-line: `--xlrg_job 200GB --lrg_job 100GB`.

### Queue size 

The queue size parameter `--queue_size 0`, is used to limit the total number of DRAM2 jobs which can occur at the same time. Thus, if the user desires to use no more than 40 CPUs at a given time, and they requested 10 CPUs for each process, the would use the command-line options: `--threads 10 --queue_size 4`. 

A queue size of 0, `--queue_size 0`, means unlimited and is the default.

------

## Example command-line options implemented in DRAM2

**Example 1: The user has unlimited computational resources and is not concerned about having to share resources BUT is using SLURM:**

`--queue_size 0 --max_forks_single_cpu 100 --max_forks_user_cpu 100 --max_forks_hmm 100`

Result: SLURM will submit as many jobs as it can as long as there are resources available.

**Example 2: The user has moderate computational resources and is concerned about having to share resources and is using SLURM:**

`--queue_size 4 --cpus 10`

Leaving the maxForks parameters alone. 

Result: SLURM will submit up to 4 jobs at a time and thus the user will, at most, be using 40 CPUs.

**Example 2: The user is only using a laptop with 16 CPUs and 16GB RAM:**

`--queue_size 2 --cpus 5`

Leaving the maxForks parameters alone. 

Instead of typing all of the memory options out on the command-line, the user modifies the `nextflow.config` to reduce the memory allotment to each process:

```bash
    /* Memory allocation */
        xsml_job = '200MB'
        sml_job = '5GB'
        med_job = '5GB'
        lrg_job = '5GB'
        xlrg_job = '5GB'
```

Result: SLURM will submit up to 2 jobs at a time and thus the user will, at most, be using 10 CPUs. Additionally, the user will at most use 10GB of RAM. This will ensure the user does not crash their laptop by saving some remaining resources for normal laptop operation.



------


## Suggestions for future development

### Memory presets

It is desired to have memory presets for use on computers with different resources. This will allow the user to tweak less parameters by hand. For example there may be a command-line option for `--memory` where there are 3 options: `laptop`, 'sml_server` and 'lrg_server'. All this would entail is some additional parameters within the `nextflow.config`.

The idea above could also be expanded to include queue size and CPU allocation.

### Profiles

It is suggested to make server-specific profiles. For example, create a W2 profile (default Singularity + SLURM) which has the partition and time (set to max) presets. This is suggested because SLURM is setup differently on most HPCs and therefore it is hard to generalize the configuration files across multiple servers. 