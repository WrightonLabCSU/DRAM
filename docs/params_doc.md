# WrightonLabCSU/dram pipeline parameters

DRAM (Distilled and Refined Annotation of Metabolism) is a tool for annotating metagenomic assembled genomes

## Pipeline Steps

Which steps to run in the pipeline.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `rename` | Rename FASTA headers based on file name. Example: sample1.fa --> (fasta header renamed to) > sample1......  Why? DRAM output is focused on scaffolds/contigs with respect to each provided input input_fasta. Thus, without renaming FASTA headers, the individual scaffolds/contigs will not be distinguashable. *If you have already renamed your FASTA headers, do not include with '--call'. | `boolean` |  |  |  |
| `call` | Whether to call genes on the input FASTA files using Prodigal. | `boolean` |  |  |  |
| `annotate` | Annotate called genes using downloaded databases. | `boolean` |  |  |  |
| `merge_annotations` | Path to directory pointing to DRAM annotation to merge into one file. This is ran as a separate pipeline. | `string` |  |  |  |
| `distill_topic` | Topic to distill. Options:  <carbon|energy|misc|nitrogen|transport> OR <default = carbon,energy,misc,nitrogen,transport> If more than one topic included, they must be seperated by a comma (--distill_topic carbon,energy or --distill_topic 'carbon,energy'). | `string` |  |  |  |
| `distill_ecosystem` | Ecosystem to distill. Options: <eng_sys|ag> If more than one ecosystem included, they must be seperated by a comma (--distill_ecosystem eng_sys,ag or --distill_ecosystem 'eng_sys,ag'). | `string` |  |  |  |
| `distill_custom` | <path/to/custom_distillate.tsv>  Custom distill file to use. | `string` |  |  |  |
| `product` | Generate a product visualization of the annotations and save the output to the output directory. | `boolean` |  |  |  |
| `format_kegg` | Format KEGG database for use in DRAM. Standalone operation, will exit after completion. | `boolean` |  |  |  |

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input_fasta` | Path to FASTA directory <details><summary>Help</summary><small>This parameter is *mandatory*.</small></details>| `string` |  |  |  |
| `fasta_fmt` | Input format for the FASTA file. | `string` | *.f* |  |  |
| `input_genes` | Directory containing called genes. Only used when not running call. This allows you to provide pre-called genes to the pipeline. | `string` |  |  |  |
| `genes_fmt` | Input format for the Genes file. | `string` | *.faa |  |  |
| `genes_fna_fmt` | Input format for the Genes fna file. Only needed if you are not running call and you are passing `--generate_gff` or `--generate_gbk`. | `string` | *.fna |  |  |

## Call Prodigal Options

Call genes on the input FASTA files using Prodigal.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `prodigal_mode` | Mode for Prodigal gene calling. | `string` | single |  |  |
| `prodigal_trans_table` | Translation table to use. | `number` | 11.0 |  |  |
| `min_contig_len` | Minimum contig length in base pairs. | `number` | 2500.0 |  |  |

## Annotate Options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `use_camper` | Use the CAMPer database for annotation. | `boolean` |  |  |  |
| `use_canthyd` | Use the Cant_Hyd database for annotation. | `boolean` |  |  |  |
| `use_dbcan` | Use the DBCan database for annotation. | `boolean` |  |  |  |
| `use_fegenie` | Use the FeGenie database for annotation. | `boolean` |  |  |  |
| `use_kegg` | Use the KEGG database for annotation. | `boolean` |  |  |  |
| `use_kofam` | Use the Kofam database for annotation. | `boolean` |  |  |  |
| `use_merops` | Use the MEROPS database for annotation. | `boolean` |  |  |  |
| `use_methyl` | Use the Methyl database for annotation. | `boolean` |  |  |  |
| `use_pfam` | Use the PFAM database for annotation. | `boolean` |  |  |  |
| `use_sulfur` | Use the Sulfur database for annotation. | `boolean` |  |  |  |
| `use_uniref` | Use the UniRef database for annotation. | `boolean` |  |  |  |
| `add_annotations` | Path to old annotations to add to the current run. | `string` |  |  |  |
| `bit_score_threshold` | Minimum BitScore of search to retain hits') | `number` | 60.0 |  |  |
| `rbh_bit_score_threshold` | Minimum BitScore of reverse best hits to retain hits | `number` | 350.0 |  |  |
| `database_list` | Database list for generating GFF and GBK. Comma sepeterated list of databases to include in the annotates. Use empty for all. Example: 'kegg,dbcan,kofam,merops,viral,camper,cant_hyd,fegenie,sulfur,methyl,uniref,pfam,vogdb' | `string` | empty |  |  |
| `generate_gff` | Generate GFF output for each input_fasta. | `boolean` |  |  |  |
| `generate_gbk` | Generate GBK output for each input_fasta. | `boolean` |  |  |  |

## Distill Options

The purpose of DRAM distill is to distill down annotations based on curated distillation summary form(s). It can be ran with either --distill_topic, --distill_ecosystem, or --distill_custom (or some combination). User's may also provide a custom distillate via --distill_custom <path/to/file> (TSV forms). Distill can be ran independent of --call and --annotate however, annotations must be provided (--annotations <path/to/annotations.tsv>). Optional tRNA, rRNA and bin quality may also be provided.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `annotations` | Path to annotations tsv to distill. Required if running without annotate | `string` |  |  |  |
| `distill_carbon_sheet` | File path to distill carbon sheet. | `string` |  |  | True |
| `distill_energy_sheet` | File path to distill energy sheet. | `string` |  |  | True |
| `distill_misc_sheet` | File path to distill misc sheet. | `string` |  |  | True |
| `distill_nitrogen_sheet` | File path to distill nitrogen sheet. | `string` |  |  | True |
| `distill_transport_sheet` | File path to distill transport sheet. | `string` |  |  | True |
| `distill_camper_sheet` | File path to distill camper sheet. | `string` |  |  | True |
| `distill_eng_sys_sheet` | File path to distill energy ecosystem sheet. | `string` |  |  | True |
| `distill_ag_sheet` | File path to distill agriculture ecosystem sheet. | `string` |  |  | True |
| `distill_dummy_sheet` | File path to distill dummy sheet. | `string` |  |  | True |

## Product Options

The purpose of DRAM --product is to generate a product visualization of the annotations and save the output to the output directory.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `groupby_column` | Column to to group by in the annotations file for etc and function groupings. | `string` | input_fasta |  |  |

## RNA Options

rRNA and tRNA input sheets, used when running DRAM distill without --call

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `rrnas` | Path to rRNA tsv. | `string` |  |  |  |
| `trnas` | Path to tRNA tsv. | `string` |  |  |  |

## Bin Quality and Taxonomy Options

Paths to bin quality and taxonomy tsvs, used after annotate and before distill

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `bin_quality` | Path to bin quality tsv. CheckM and CheckM2 compatible. | `string` |  |  |  |
| `taxa` | Path to bin taxonomy tsv. Compatible with GTDB. | `string` |  |  |  |

## Database Options

File paths to databases used in the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `kegg_db` |  | `string` | ${launchDir}/databases/kegg/ |  | True |
| `uniref_db` |  | `string` | ${launchDir}/databases/uniref/ |  | True |
| `pfam_mmseq_db` |  | `string` | ${launchDir}/databases/pfam/mmseqs/ |  | True |
| `merops_db` |  | `string` | ${launchDir}/databases/merops/ |  | True |
| `viral_db` |  | `string` | ${launchDir}/databases/viral/ |  | True |
| `kofam_db` |  | `string` | ${launchDir}/databases/kofam/ |  | True |
| `kofam_list` |  | `string` | ${launchDir}/databases/kofam/kofam_ko_list.tsv |  | True |
| `dbcan_db` |  | `string` | ${launchDir}/databases/dbcan/ |  | True |
| `dbcan_fam_activities` |  | `string` | ${launchDir}/databases/dbcan/dbcan.fam-activities.tsv |  | True |
| `dbcan_subfam_activities` |  | `string` | ${launchDir}/databases/dbcan/dbcan.fam-activities.tsv |  | True |
| `vog_db` |  | `string` | ${launchDir}/databases/vog/ |  | True |
| `vog_list` |  | `string` | ${launchDir}/databases/vogdb/vog_annotations_latest.tsv.gz |  | True |
| `camper_hmm_db` |  | `string` | ${launchDir}/databases/camper/hmm/ |  | True |
| `camper_hmm_list` |  | `string` | ${launchDir}/databases/camper/hmm/camper_hmm_scores.tsv |  | True |
| `camper_mmseqs_db` |  | `string` | ${launchDir}/databases/camper/mmseqs/ |  | True |
| `camper_mmseqs_list` |  | `string` | ${launchDir}/databases/camper/mmseqs/camper_scores.tsv |  | True |
| `canthyd_hmm_db` |  | `string` | ${launchDir}/databases/canthyd/hmm/ |  | True |
| `cant_hyd_hmm_list` |  | `string` | ${launchDir}/databases/canthyd/hmm/cant_hyd_HMM_scores.tsv |  | True |
| `canthyd_mmseqs_db` |  | `string` | ${launchDir}/databases/canthyd/mmseqs/ |  | True |
| `canthyd_mmseqs_list` |  | `string` | ${launchDir}/databases/canthyd/mmseqs/cant_hyd_BLAST_scores.tsv |  | True |
| `fegenie_db` |  | `string` | ${launchDir}/databases/fegenie/ |  | True |
| `fegenie_list` |  | `string` | ${launchDir}/databases/fegenie/fegenie_iron_cut_offs.txt |  | True |
| `sulfur_db` |  | `string` | ${launchDir}/databases/sulfur/ |  | True |
| `methyl_db` |  | `string` | ${launchDir}/databases/methyl/ |  | True |
| `sql_descriptions_db` |  | `string` | ${launchDir}/databases/db_descriptions/description_db.sqlite |  | True |
| `kegg_e_value` |  | `string` | 1e-05 |  | True |
| `kofam_e_value` |  | `string` | 1e-05 |  | True |
| `dbcan_e_value` |  | `string` | 1e-15 |  | True |
| `merops_e_value` |  | `string` | 1e-1 |  | True |
| `vog_e_value` |  | `string` | 1e-05 |  | True |
| `camper_e_value` |  | `string` | 1e-05 |  | True |
| `uniref_e_value` |  | `string` | 1e-05 |  | True |
| `canthyd_e_value` |  | `string` | 1e-05 |  | True |
| `sulfur_e_value` |  | `string` | 1e-05 |  | True |
| `fegenie_e_value` |  | `string` | 1e-05 |  | True |

## Format Kegg Options

Options for preparing the KEGG database for use in DRAM.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `kegg_pep_root_dir` | Only required if you need to concatenate all of KEGG's provided pep files. Root directory to downloaded KEGG peptide files. The pipeline will search for all pep files in this directory and concatenate them into a single file. In format <root_dir>/*/*.pep. | `string` |  |  |  |
| `kegg_pep_loc` | Path to and of the gene fasta files that are provided by the KEGG FTP server or a concatenated version of them. Either this or kegg_pep_root_dir must be provided. | `string` |  |  |  |
| `gene_ko_link_loc` | Path to and of the KO file that is provided by the KEGG FTP server. | `string` |  |  |  |
| `skip_gene_ko_link` | Skip the gene_ko_link file. If you are using an older version of KEGG that does not supply the gene_ko_link file you can use this option to skip the gene_ko_link file. Otherwise, the gene_ko_link file is required. | `boolean` |  |  |  |
| `kegg_download_date` | The date the KEGG database was downloaded. If not provided, the current date will be used. | `string` | yyyy-MM-dd |  |  |

## SLURM Options

Generic options for SLURM job submission. More customized options can be configured in your own Nextflow config file. See https://www.nextflow.io/docs/latest/executor.html#slurm and example configs here: https://github.com/nf-core/configs/tree/master/conf

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `slurm` | Launch the pipeline using the SLURM executor. Without this option, the pipeline will run in your current shell/environment. | `boolean` |  |  |  |
| `partition` | Name of the SLURM partition to use for job submission. If not provided, the default partition will be used. | `string` |  |  |  |
| `queue_size` | Maximum number of jobs to submit to the queue at once. | `integer` | 10 |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `threads` |  | `integer` | 10 |  | True |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  | True |
| `trace_report_suffix` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. | `string` | yyyy-MM-dd_HH-mm-ssZ |  | True |
