Resource Management
===================

DRAM2 is implemented in Nextflow.

Nextflow is able to scale horizontally on a given machine.

What does this mean?

Horizontal scaling refers to the ability to distribute computational tasks across multiple computing resources, such as cores or nodes, in parallel. In the context of Nextflow, this means that a single workflow can leverage the computational power of multiple CPUs or nodes, allowing for faster execution of tasks and improved overall performance.

By utilizing horizontal scaling, Nextflow can efficiently manage and execute workflows that require significant computational resources, such as those involved in genomic data analysis. This enables DRAM2 to process large datasets and complex analyses in a timely manner, making it suitable for a wide range of research and bioinformatics applications.

Summary
--------

DRAM2 comes with configuration files which have the option to change how many "things" can happen at a time in the pipeline.

A user can modify these, "maxForks", parameters within the ``nextflow.config`` to increase the number of "things" which DRAM2 can perform at a given time.

**NOTE**: Development is in progress to enable different DRAM2 modes: "lite", "medium" and "heavy". Where "lite" would be for a good laptop and "heavy" for a HPC. These options will alter the CPU and memory (RAM) requirements for each process.

---------------

CPU Management
---------------

CPU management in Nextflow can be tricky BUT it gives the user great power to scale up their DRAM2 run.

Why is it tricky?

Each process in the DRAM2 pipeline is submitted as a separate job.

Thus, each separate annotation for each sample, for each database, is submitted as a separate job.


Nextflow maxForks parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Nextflow, "maxForks" refers to the maximum number of concurrent tasks that can be executed simultaneously. It is a parameter that you can specify when running a Nextflow workflow to control the level of parallelism.\

Imagine you have a workflow that consists of multiple tasks, each performing a different computation or analysis. By default, Nextflow will try to execute as many tasks as possible concurrently, depending on the available computational resources (such as CPU cores or nodes). However, setting a "maxForks" value allows you to limit the maximum number of tasks that can run simultaneously.

For example, if you set "maxForks" to 4, Nextflow will only allow up to 4 tasks to run concurrently at any given time, even if there are more tasks ready to be executed. This can be useful for controlling resource usage and preventing overloading of the system, especially in environments with limited computational resources.

In simple terms, "maxForks" helps you manage the level of parallelism in your workflow by controlling how many tasks can be executed at the same time. It allows you to balance the trade-off between maximizing throughput and avoiding resource contention.

DRAM2 has 3 types of processes:

1. small processes which can only use 1 CPU
2. HMM searching with can only use 2 CPUs
3. All other processes which can utilize more than 2 CPUs

The "maxForks" for each type of process is set in the ``nextflow.config``:

    .. code-block:: bash

        /* Max forks Options */
            max_forks_single_cpu  = 2
            max_forks_user_cpu = 2
            max_forks_hmm = 2

---------------

Memory Management
------------------

**Memory management is only utilized when the SLURM-based profile option is used: ``-profile conda_slurm`` or ``-profile singularity_slurm``.**

By default, DRAM2 has memory parameters set at what we feel is required to run DRAM2 at its best on a small HPC with 512GB RAM and 96 CPUS.

These are as follows within the ``nextflow.config``:

    .. code-block:: bash

        /* Memory allocation */
            xsml_job = '200MB'
            sml_job = '12GB'
            med_job = '50GB'
            lrg_job = '200GB'


**Note:** These are minimum memory allocations for SLURM and can be modified by the user either in the nextflow.config file or by using command-line options in a DRAM2 run: ``--lrg_job 100GB``.

---------------

Examples
--------------- 

Example 1:
~~~~~~~~~~

John Doe is computing on a MAC laptop w/ 8-core CPU and 16GB RAM

maxForks: John modifies the ``nextflow.config`` and sets all maxForks = 1:

    .. code-block:: bash

        /* Max forks Options */
            max_forks_single_cpu  = 1
            max_forks_user_cpu = 1   (provided by --threads <number>)
            max_forks_hmm = 1

CPUs: John wants to make sure he can still use his computer so he only uses 6 CPUs: ``--threads 6``.

While DRAM2 will want to submit more jobs at a time, the additional jobs will have to wait until the previous job has stopped.

John also excluded the following databases because of their size and memory requirements: KEGG and UniRef.



Example 2:
~~~~~~~~~~

Jane Doe is computing on a HPC with 20 nodes (each node has 96 CPUs and 512GB RAM). The HPC also requires the use of SLURM.

maxForks: Jane modifies the ``nextflow.config`` and allows 5 single-CPU jobs to run at a time, allows 10 jobs to run at a time on the ``--threads`` she specifies, and allows 20 HMM searches to occur at a time.

    .. code-block:: bash

        /* Max forks Options */
            max_forks_single_cpu  = 5
            max_forks_user_cpu = 10   (provided by --threads <number>)
            max_forks_hmm = 50

CPUs: Jane wants to make sure DRAM2 runs efficiently across the HPC: ``--threads 40``.

DRAM2 will submit as many jobs as there are resources. Remaining jobs will be places in the SLURM queue and will be performed when resources are available.

Jane included all of the databases to annotate from.

Outcome:

Say Jane had 1000 samples. This means, because of the maxForks of ``max_forks_user_cpu = 10``, 10 samples will be annotated at the same time for the KEGG database, each with 40 CPUs = 400 CPUs.

Jane was a bit ambitious - this HPC cannot possibly run all annotations, for all databases, for all samples, at the same time. However, because this HPC uses SLURM, her jobs will sit in the queue and will be ran when resources are available. 