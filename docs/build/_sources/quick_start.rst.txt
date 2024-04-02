Quick Start Guide
=================

Assuming you have followed the :doc:`installation` instructions and you have the databases you desire, you can start running the commands below:

DRAM2 apps Call, Annotate, and Distill can all be run at once, or alternatively, each app can be run individually (assuming you provide the required input data for each app..)

Additionally, `--merge-annotations` and `--rename` can be run independently of any other apps.

Example command-line usage
---------------------------

1. Rename fasta headers based on input sample file names:

    .. code-block:: bash

        nextflow run DRAM2.nf --rename --input_fasta_dir <path/to/fasta/directory/>

2. Call genes using input fastas (use --rename to rename FASTA headers.:

    .. code-block:: bash

        nextflow run DRAM2.nf --call --rename --input_fasta_dir <path/to/fasta/directory/>

3. Annotate called genes using input called genes and the KOFAM database:

    .. code-block:: bash

        nextflow run DRAM2.nf --annotate --input_genes <path/to/called/genes/directory> --use_kofam

4. Annotate called genes using input fasta files and the KOFAM database:

    .. code-block:: bash

        nextflow run DRAM2.nf --annotate --input_fasta <path/to/called/genes/directory> --use_kofam

5. Merge various existing annotations files together (Must be generated using DRAM2..:

    .. code-block:: bash

        nextflow run DRAM2.nf --merge_annotations <path/to/directory/with/multiple/annotation/TSV/files>

6. Distill using input annotations:

    .. code-block:: bash

        nextflow run DRAM2.nf --distill_<topic|ecosystem|custom> --annotations <path/to/annotations.tsv>

7. (Combined): Call, annotate and distill input fasta files:

    .. code-block:: bash

        nextflow run DRAM2.nf --rename --call --annotate --use_<database(s) --distill_topic <distillate(s.>

8. Call and Annotate genes using input fastas and KOFAM database. Distill using the default topic and the AG ecosystem:

    .. code-block:: bash

        nextflow run DRAM2.nf --input_fasta_dir <path/to/fasta/directory/> --outdir <path/to/output/directory/> --call --annotate --distill_topic default --distill_ecosystem ag --threads <threads> --use_kofam

9. "Real-world" example using the test data provided in this repository:

    .. code-block:: bash

        nextflow run -bg DRAM2.nf --input_fasta ../test_data/DRAM2_test_data/ --outdir DRAM2-test-data-call-annotate-distill --threads 8 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom test-data/custom-test-distilalte.tsv --profile conda_slurm --slurm_node main -with-report -with-trace -with-timeline

Breakdown of example (9):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- `-bg` : Nextflow option to push the run immediately into the background. (Thus, you can log out on an HPC and the run will continue..
- `-profile` : Nextflow option to select profile (Conda vs Singularity and SLURM vs no-SLURM..
- `--slurm_node`: DRAM2 option to select a specific node to compute on during the whole run.
- `-with-trace`: Nextflow option to output a process-by-process report of the run. (TEXT.
- `-with-report`: Nextflow option to output a process-by-process report of the run. (HTML.
- `-with-timeline`: Nextflow option to output a process-by-process HTML timeline report of the run. (HTML.


Select a Profile which suits your compute setup
------------------------------------------------

DRAM2 utilizes either Conda or Singularity for dependency management and the user MUST choose one of the following options on execution of any DRAM2 command

*The Nextflow profile option is used (`-profile`) - yes! a single hyphen!*

1. 

    .. code-block:: bash
    
        -profile conda


    This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

2. 

    .. code-block:: bash
    
        -profile conda_slurm

    This option will submit each individual DRAM2 process as its own SLURM job. (See Wiki Resource Management for details).
    This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

3. 

    .. code-block:: bash
    
        -profile singularity

    This option relies on the local systems Singularity. Nextflow will create its own Conda environments to run in. 

4. 

    .. code-block:: bash
    
        -profile singularity_slurm

    This option will submit each individual DRAM2 process as its own SLURM job. (See Wiki Resource Management for details)
    This option relies on the local systems Singularity to run the downloaded Singularity container.  