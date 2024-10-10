# DRAM2 Nextflow Profiles

To enable flexible computing across different computers/HPCs/Servers/etc, Nextflow provides the ability to create custom profiles which users can specify.

Explanation by example:

Some HPCs, like W2, Riviera and Alpine, require users to submit jobs using SLURM. However, some DRAM2 users may want to run DRAM2 on their laptop and have no need for SLURM scheduling. Thus, DRAM2 provides different profiles for these various computing types.


## How profiles are chosen by the user

*The Nextflow profile option is used (`-profile`) - yes! a single hyphen!*

1) `-profile conda`
   
  This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

2) `-profile conda_slurm`
   
  This option will submit each individual DRAM2 process as its own SLURM job. (See Wiki Resource Management for details).
  This option relies on the local systems Conda. Nextflow will create its own Conda environments to run in. 

3) `-profile singularity`
   
  This option relies on the local systems Singularity. This option relies on the local systems Singularity to run the downloaded Singularity container.  

4) `-profile singularity_slurm`
   
  This option will submit each individual DRAM2 process as its own SLURM job. (See Wiki Resource Management for details)
  This option relies on the local systems Singularity to run the downloaded Singularity container.  
