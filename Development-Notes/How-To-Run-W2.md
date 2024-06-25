
W2 location to run DRAM2:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF`

Super small - Call, Annotate, Distill:
Has -resume for illustration:
```groovy
nextflow run DRAM2.nf --input_fasta /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/subsampled/super_small/ --outdir DRAM2-super-small-test-2-06252024 --threads 10 --call --rename --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees -resume
```


Super small - Call, Annotate, Distill:
NO resume:
```groovy
nextflow run DRAM2.nf --input_fasta /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/subsampled/super_small/ --outdir DRAM2-super-small-test-2-06252024 --threads 10 --call --rename --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees
```


Super small - Starting from Annotate:
```groovy
nextflow run DRAM2.nf --input_genes DRAM2-super-small-test-2-06252024/Prodigal_v2.6.3/ --outdir DRAM2-super-small-test-2-ANNOTATE-06252024 --threads 10 --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees 
```


Super small - Starting from Distill:
```groovy
nextflow run DRAM2.nf --annotations DRAM2-super-small-test-2-06252024/RAW/raw-annotations.tsv --outdir DRAM2-super-small-test-2-DISTILL-06252024 --threads 10 --annotate --use_camper --use_canthyd --distill_topic default --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --no_trees 
```


Bigger, subsampled - Call, Annotate, Distill:
```groovy
nextflow run DRAM2.nf --input_fasta /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/subsampled/ --outdir DRAM2-super-small-test-2-06252024 --threads 5 --call --rename --annotate --use_kofam --use_kegg --use_dbcan --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 0 -profile singularity_slurm -resume --no_trees
```


