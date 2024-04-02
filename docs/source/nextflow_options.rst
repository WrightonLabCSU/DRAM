Nextflow Pipeline Options
=========================

The DRAM2 pipeline offers a range of options that leverage Nextflow's capabilities to enhance the flexibility and efficiency of your workflow. This page details some of the core Nextflow options you can use when running the DRAM2 pipeline.

Example Command
---------------

.. code-block:: bash

    nextflow run DRAM2.nf --rename --call --annotate --use_<database(s)> --distill_topic <distillate(s)>

Options Breakdown
^^^^^^^^^^^^^^^^^

1. 

    .. code-block:: bash
    
        -bg

    Pushes the run immediately into the background. This is particularly useful on HPC environments, allowing you to log out while the run continues.

2. 

    .. code-block:: bash
    
        -profile <profile>

    Allows the selection of a specific profile (e.g., Conda vs Singularity, SLURM vs no-SLURM). Profiles help manage different execution environments and resource management systems.

3. 

    .. code-block:: bash
    
        -with-trace

    Outputs a detailed process-by-process report of the run in TEXT format.

4. 

    .. code-block:: bash
    
        -with-report

    Generates a comprehensive process-by-process report of the run in HTML format.

5. 

    .. code-block:: bash
    
        -with-timeline

    Produces a process-by-process HTML timeline report of the run, offering a visual representation of the pipeline's execution over time.

Nextflow Tips and Tricks
------------------------

The `-resume` option
^^^^^^^^^^^^^^^^^^^^^^^^
In Nextflow DSL2, the ``-resume`` option provides a powerful mechanism for efficiently managing workflow runs. It allows you to resume a run, modify parameters, and reuse previously generated data. Here are common scenarios where ``-resume`` is beneficial:


Scenario: Adding Additional Databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you want to run DRAM2 again but with additional annotation databases. Using the ``-resume`` option enables you to reuse your called genes and existing annotations, provided the ``work/`` directory has not been deleted.

For instance, if your initial run used ``--use_kofam --use_dbcan``, you can add more databases like ``--use_kegg --use_uniref`` to your command. This approach reuses the existing called genes to annotate with the newly added databases while retaining the existing annotations, assuming the modifications are retained in the command.

Limitations of the `-resume` Option
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It's important to note how the ``-resume`` option does not work:

- The option does not allow for rerunning a workflow from scratch without manually deleting the ``work/`` directory or using a different project directory. It is designed to pick up where a previous run left off, utilizing cached results to avoid redoing work.

This documentation aims to clarify the usage and benefits of Nextflow's options within the DRAM2 pipeline, enhancing your workflow's efficiency and flexibility.
