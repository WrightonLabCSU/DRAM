
# Reed Woyda daily logs

Note:
This may or may not be useful for a future developer. This document contains most daily happenings and development notes for DRAM2. If anything, it contains an expansive list of test commands. 

## Friday October 6

- Need to test kegg-kofam database
	- Still have to write this into the main nf script
	- finish writing mmseqs2.nf and incorporate into the main workflow
	- write 
- I may need to re-build the kegg (and other) databases to ensure I know what format they are in. 
	- Ran into this error:
```
Caused by:
  Process `MMSEQS2 (BB02)` terminated with an error exit status (1)

Command executed:

  mkdir tmp
  
  
  # make fasta query into database
  # Creates queryDB, queryDB_h, queryDB.index, queryDB_h.index and queryDB.lookup
  mmseqs createdb BB02_called_genes.fna queryDB
  
  # make query to target db
  
  mmseqs search queryDB kegg.20220928.mmsdb query_target_db tmp --threads 20

Command exit status:
  1

Command output:
  Mask lower case residues               	0
  Minimum diagonal score                 	15
  Spaced k-mers                          	1
  Spaced k-mer pattern                   	
  Local temporary path                   	
  Rescore mode                           	0
  Remove hits by seq. id. and coverage   	false
  Sort results                           	0
  Mask profile                           	1
  Profile E-value threshold              	0.1
  Global sequence weighting              	false
  Allow deletions                        	false
  Filter MSA                             	1
  Maximum seq. id. threshold             	0.9
  Minimum seq. id.                       	0
  Minimum score per column               	-20
  Minimum coverage                       	0
  Select N most diverse seqs             	1000
  Min codons in orf                      	30
  Max codons in length                   	32734
  Max orf gaps                           	2147483647
  Contig start mode                      	2
  Contig end mode                        	2
  Orf start mode                         	1
  Forward frames                         	1,2,3
  Reverse frames                         	1,2,3
  Translation table                      	1
  Translate orf                          	0
  Use all table starts                   	false
  Offset of numeric ids                  	0
  Create lookup                          	0
  Add orf stop                           	false
  Overlap between sequences              	0
  Sequence split mode                    	1
  Header split mode                      	0
  Chain overlapping alignments           	0
  Merge query                            	1
  Search type                            	0
  Search iterations                      	1
  Start sensitivity                      	4
  Search steps                           	1
  Exhaustive search mode                 	false
  Filter results during exhaustive search	0
  Strand selection                       	1
  LCA search mode                        	false
  Disk space limit                       	0
  MPI runner                             	
  Force restart with latest tmp          	false
  Remove temporary files                 	false

Command error:
  Input database "kegg.20220928.mmsdb" has the wrong type (Generic)
  Allowed input:
  - Index
  - Nucleotide
  - Profile
  - Aminoacid

Work dir:
  /home/rwoyda/Projects/DRAM2-Nextflow/DRAM2-NF/work/a8/42fb4c31870876bc63e431901dd2d2


```


## Tuesday October 17
- KEGG is proving quite tricky
	- I have tried
		- Loading in the pre-built KEGG mmsdb used from DRAM
			- Keep getting the error:
```
Command error:
  Input database "index/index" has the wrong type (Generic)
  Allowed input:
  - Index
  - Nucleotide
  - Profile
  - Aminoacid

```

- Also have tried to build an index from this but that did not work either
- Also tried to pass in already build indices but still could not search them

For now - will come back to the mmseqsdb

Going to switch it up and go to HMMER. 
- In the `DRAM2-Nextflow-Database-recipe-v3` (DRAM2-NF-Databases-kofam.sif)
	- We have `kofam_profiles.hmm*` files.
	- Going to write the hmm search processes and try to annotate this way.
	- I wrote multiple things
		- ` kegg_hmm_formatter.nf`
			- `depends on kegg_hmm_formatter.py`
			- Tested this with kofam - works with kofam!
		- generic_hmm_formatter.nf
			- `generic_hmm_formatter.py`
			- Not tested
			- Want to bring in another database to test this one
		- looks like a lot of the databases have custom hmm_formatters
			- dbcan has `dbcan_hmmscan_formater`
			- vogdb has own as well: `vogdb_hmmscan_formater`
		- Do these all have to be different!?!?!
			- Okay if it is, just makes so many different functions which are very similar..
	- Created this image for further testing
		- `DRAM2-NF-Databases-kofam-vogdb-dbcan.sif`
			- Testing now - WORKS
	- Wrote distill processes
		- need to test - NOT DONE - but going to change anyway to deal with combined_annotations output
	- Wrote combine annotations process
		- first iteration of the python script to merge the various annotations files
			- need to test - NOT DONE
			- need to generalize the inputs - i did some but may need more
		- need to add in the two columns: sample and contig/scaffold - where to get this info from? from prodigal?
			- may need to combine the prodigal output and parse it to add these values
		- need to test



## Wednesday October 18
- Fix rename fasta
	- This relies on bbtools
	- For now, just going to use a Gene-Centric container -DONE
	- Need to update the container - NOT DONE
- Today I am going to focus on getting combine_annotations and the distill steps working
	- The current distill needs to be changed to deal with the new combined_annotations output
- Ideal next steps would be to
	- Get an mmseqdb working - maybe merops?
	- Get heatmap HTML and tsv working



### Rebuild singularity image
- `DRAM2-NF-Databases-kofam-vogdb-dbcan.sif`
	- need to include all dbcan database files - I needed a * after the path/to/dbCAN-HMMdb-V11.txt
	- include bbtools
	- include dependencies for the python scripts
		- pandas - already in here it looks like!
	- REBULDING NOW
		- `DRAM2-NF-Databases-kofam-vogdb-dbcan-v4.sif DRAM2-Nextflow-Database-recipe-v4`
		- NOT TESTED
		- Need to add openpyxl
- Distill summary works!
- Need to work on distill final and then the product heatmap
	- Distill final works! Even with multiple modules like engineered systems!!
	- Still need to do heatmap and test this on samples with different names
		- started working on the heatmap but the truth is I have no idea what is going on - Need to ask Kayla how this is performed
		- `nextflow run DRAM2.nf --input_fasta_dir test_data/subsampled/small/ --outdir DRAM2-kofam-distill_summary-10182023 --call --rename --annotate --use_kofam --distill --add_annotations Engineered_systems --threads 20 -resume`
- Going to switch to getting another database to work - will try dbcan
- 
```
great, we need to use two additional input files (TSV) to add additional columns to the output.
We want to include these columns in the output:
dbcan_ids	dbcan_hits	dbcan_subfam_ec	dbcan_best_hit


First new input file:
${params.dbcan_fam_activities}
Looks like this:
# Copyright disclosure: This data is parsed from webpages of www.cazy.org.	
AA0	  
AA10	  AA10 (formerly CBM33) proteins are copper-dependent lytic polysaccharide monooxygenases (LPMOs); some proteins have been shown to act on chitin, others on cellulose; lytic cellulose monooxygenase (C1-hydroxylating) (EC 1.14.99.54); lytic cellulose monooxygenase (C4-dehydrogenating)(EC 1.14.99.56); lytic chitin monooxygenase (EC 1.14.99.53); lytic xylan monooxygenase / xylan oxidase (glycosidic bond-cleaving) (EC 1.14.99.-)
AA11	  AA11 proteins are copper-dependent lytic polysaccharide monooxygenases (LPMOs); cleavage of chitin chains with oxidation of C-1 has been demonstrated for a AA11 LPMO from Aspergillus oryzae;lytic chitin monooxygenase (1.14.99.53)
AA12	  The pyrroloquinoline quinone-dependent oxidoreductase activity was demonstrated for the CC1G_09525 protein of Coprinopsis cinerea. 
AA13	  AA13 proteins are copper-dependent lytic polysaccharide monooxygenases (LPMOs); cleavage of starch with oxidation of C-1 at the site of cleavage has been demonstrated for the LPMO encoded by gene NCU08746 from Neurospora crassa; lytic starch monooxygenase / starch oxidase (glycosidic bond-cleaving) (EC 1.14.99.55)

It has two columns. The first column contains the "dbcan_ids".

The second input file is:
```





## Tuesday November 21
- Looks like this will be a real thing! 
	- Starting in January
- Want pre-build databases NOT within containers
	- Long term - have a repo online where people can download the pre-built databases (NOT KEGG)



## Thursday November 30
- Worked on organizing, cleaning and documenting DRAM2.nf code
- Building new Singularity container w/o databases and forms
	- `DRAM2-Nextflow-No-Databases-Nov302023-V1`
	- Put forms in `assets/forms/`
- Testing call: `nextflow run DRAM2.nf --input_fasta test_data/subsampled/super_small/ --outdir DRAM2-supersmall-Call-test1 --threads 10 --call --rename`
	- WORKED!
- Mybe test annotate tomorrow with kofam


## Thursday December 28
- Boy has it been awhile..
- Trying to get my head wrapped around where DRAM2 in nextflow sits.
- To recap:
	- DRAM2 NF code was re-written/structured to set a good foundation for the furture
	- 

## Wednesday January 3
- Pulled DRAM2 NF GitHub to Riviera
- Did simple call test:
```groovy
nextflow run DRAM2.nf --input_fasta test_data/subsampled/super_small/ --outdir DRAM2-supersmall-Call-test1 --threads 50 --call --rename
```
- Worked!
- Need to test annotate with the kofam database



## Monday January 8
- I did not do a good job last week of recording my work here
	- To recap:
		- DBCAN works!
			- Has its own dedicated formatter which uses both the fam and sub fam activities spreadsheets
			- DOES not have "top-hit" functionality yet
		- KOFAM works!
			- Using the KEGG_formatted
			- Currently changing this to has its own dedicated formatter - NEED to TEST 
		- DBCAN brought in additional columns which the current COMBINE_ANNOTATIONS does not work with
			- Need to re-work the COMBINE_ANNOTATIONS and COUNT_ANNOTATIONS to work with DBCAN
				- Also, needs to be somewhat generalized
		- All of the support python scripts are hard coded within the main script - will be put in the containers in the future
Command I was testing with:
``` groovy
nextflow run DRAM2.nf --input_fasta test_data/subsampled/super_small/ --outdir DRAM2-supersmall-Call-test1 --threads 10 --call --rename --annotate --use_kofam --use_dbcan --slurm_node zenith
```



------

## Tuesday January 16
- Re-worked the #TODO section and outlined the changes to annotations and distill
	- Also added in tRNA and rRNA searching
	- Also added in ingest quality and GTDB
- Got the test data from the DRAM1 paper
	- First, going to subset this data and run it on the current version of DRAM2-NF
		- Subsampled bin 101 and bin 115:
		- `reformat.sh in=bin.101.fa out=subsampled/bin.101.fa samplerate=0.8`
	
``` groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_samples/subsampled/ --outdir DRAM2-HMP-2-subsampled --threads 10 --call --rename --annotate --use_kofam --use_dbcan --slurm_node zenith
```

Results:
```
DRAM2 Nextflow
===================================
fastas       : ../test_data/HMP_DRAM1_test_data/2_samples/subsampled/
outdir       : DRAM2-HMP-2-subsampled
threads      : 10
call genes   : true
annotate     : true
distill      : false

executor >  slurm (18)
[ac/adba2f] process > RENAME_FASTA (bin)        [100%] 2 of 2 ✔
[61/47a83e] process > CALL_GENES (bin)          [100%] 2 of 2 ✔
[1a/46fb87] process > HMM_SEARCH_KOFAM (bin)    [100%] 2 of 2 ✔
[12/4dc6dd] process > PARSE_HMM_KOFAM (bin)     [100%] 2 of 2 ✔
[5b/29b6f1] process > KOFAM_HMM_FORMATTER (bin) [100%] 2 of 2 ✔
[bb/322ff8] process > HMM_SEARCH_DBCAN (bin)    [100%] 2 of 2 ✔
[10/346aad] process > PARSE_HMM_DBCAN (bin)     [100%] 2 of 2 ✔
[ef/f13df5] process > DBCAN_HMM_FORMATTER (bin) [100%] 2 of 2 ✔
[16/73d13f] process > COMBINE_ANNOTATIONS       [100%] 1 of 1 ✔
[0a/30a87d] process > COUNT_ANNOTATIONS         [100%] 1 of 1 ✔
Completed at: 16-Jan-2024 11:51:57
Duration    : 16m 7s
CPU hours   : 0.9
Succeeded   : 18
```

Used `-resume` and added `--distill`
Results:
- Worked!
- However the sample names got truncated to just "bin" - resulting in both samples having the same name.
	- Need to fix
	- I modified the code for setting the sampleName - NEED TO TEST

- Next is to implement tRNA and rRNA searching
	- tRNA first
		- tRNAscan-SE
		- Build into the main singularity image
		- New image will be `DRAM2-Nextflow-No-Databases-Nov302023-V2(.sif)`
			- DONE
	- tRNA and rRNA both run
		- working on formatting

------
## Wednesday January 17
- Working on the tRNAscan-SE formatting
	
``` groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/1_large/ --outdir DRAM2-HMP-2-subsampled --threads 10 --call --rename --annotate --slurm_node zenith -resume

```
- DONE! I have tested it a bit but will need more.
  - Also, want to use multithreading in the future
- Working on rRNA formatting now - DONE
- Both formatters work mostly as intended. 
	- Need to keep track of pseudos
- Working on forming the distill sheets


------

## Thursday January 18
- Ran big test for tRNA collect
	- ran on the 76 HMP samples from the DRAM1 paper
Command:
`nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-2-subsampled --threads 15 --call --rename --annotate --slurm_node zenith -resume -with-report -with-trace -with-timeline`



Only one sample ran into this error:

```
(base) zenith /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF: cat work/03/51d4eb2a19ee6125673f5c42451cb0/.command.log 

tRNAscan-SE v.2.0.12 (Nov 2022) - scan sequences for transfer RNAs
Copyright (C) 2022 Patricia Chan and Todd Lowe
                   University of California Santa Cruz
Freely distributed under the GNU General Public License (GPLv3)

------------------------------------------------------------
Sequence file(s) to search:        bin-88_renamed.fna
Search Mode:                       General
Results written to:                bin-88_trna_out.txt
Output format:                     Tabular
Searching with:                    Infernal single-pass scan
                                   Fast mode
Covariance model:                  /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm

Temporary directory:               /tmp
------------------------------------------------------------

Status: Running Infernal analysis
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_2886. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_2886
Failed to open sequence file /tmp/tscan87.trna for reading

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3009. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3009
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3070. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3070
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3080. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3080
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3095. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3095
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3281. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3281
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3317. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3317
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3909. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3909
Sequence file /tmp/tscan87.trna is empty or misformatted

Error: Infernal cmsearch cannot be completed successfully for bin-88_scaffold_3997. /opt/miniconda/bin/cmsearch -g --mid --notrunc --cpu 15 /opt/miniconda/lib/tRNAscan-SE/models/TRNAinf.cm /tmp/tscan87.trna > /tmp/tscan87_cm_general.out (Exit Code: 256 Inappropriate ioctl for device)
Error: Fail to run infernal for bin-88_scaffold_3997
```


- I re-ran the run without this sample because i jsut want to test the tRNA collect
	- This worked but, in comparing the values, they do not line up. I will need to work backwards
	- This seems somewhat random. 
		- When I reran the 75 (minus bin-88), then It threw an error on bin-1
	- I re-built the singularity container and removed the install leftover from DRAM1: tRNAscan-SE version 2.0.11 and only kept 2.0.12
- There are some discrepencies between the DRAM1 published results and by results. 
	- My results do appear to be correct. Maybe this is dues to the version change?

Running full HMP test 
- validate rRNA collect and tRNA collect

Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-test-rrna --threads 15 --call --rename --annotate -with-trace -resume
```

Results:
Did not work.

Re-did with only 10 samples:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_10/ --outdir DRAM2-HMP-test-rrna-10 --threads 15 --call --rename --annotate -with-trace -resume
```
- This worked but tRNAscan-SE is very unpredictable!
	- I started and resumed this 3 times beofre it finally completed...
		- tRNAscan-SE gets an error that files already exist and it wants into to overwrrite or not.
- I am going to proceed with other stuff and them come back to these. 
	- This may be a solution: https://github.com/UCSC-LoweLab/tRNAscan-SE/issues/2


Next test: Run kofam and dbcan with the tRNA and rRNA

Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/small/ --outdir DRAM2-small-test --threads 15 --call --rename --annotate --use_kofam --use_dbcan -with-trace -resume
```
- Keep running into the tRNAscan-SE error 
	- Will comment this out for now and then proceed
- Added in flags for Bin Quality and Taxonomy: `--bin_quality AND --taxa`
	- Added file check for these on line 324
	- Added check for these flags and placeholder code for new processes in annotate workflow at line 568


------


## Wednesday January 24
- Ran GTDB, checkm1 and checkm2 on the HMP-76 samples:
	- `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/HMP_DRAM1_test_data/`
	- `checkm2-HMP-76-quality.tsv`
	- `gtdbtk.bac120.summary.tsv`
	- `checkm1-HMP-76-quality.tsv`
- Next is to write the scripts to add the necessary columns to the annotations file
	- comment out tRNA and rRNA for now

First test for add_bin_quality:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-bin-qual-test --threads 15 --call --rename --annotate --use_dbcan -with-trace -resume
```

Worked! For both checkm and checkm2

- Now working on add_taxa
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-bin-qual-test --threads 15 --call --rename --annotate --use_dbcan -with-trace -resume
```

Worked! 

- Next is to add in the start, end and strandedness cols to annotations
	- Did for dbcan and it works!
	- Modified code for kofam but need to test - on a small dataset
		- Worked!

------

## Thursday January 25
- Going to dive into the distill
	- First, need to change how distill topics and ecosystems can be added
		- Ideas:
			- User can choose to run default distill by using the flag `--distill_topic default`
			- The use can add additional topic and/or ecosystem distill by using `--distill_topic name1 name2...` and/or `--distill_ecosystem name1 name2`
			- The user can also provide custom distillate sheets with paths using the flag `--distill_custom <path> <path>....<path>`
	- Then, get the Genome Stats sheet and the tRNA and rRNA sheets added


`--distill_topic`:
- takes in one or more of these: 


`--distill_ecosys`:



`--distill_custom:


First test:
use the custom sheet: `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/distill_ecosystems/test-2-topic.tsv`
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-bin-qual-test --threads 15 --call --rename --annotate --use_kofam -with-trace --distill_topic default energy misc nitrogen transport --distill_ecosys eng_sys ag --distill_custom /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/distill_ecosystems/test-2-topic.tsv -resume 
```

-----

## Friday January 26
- Re-worked distill to use the multiple input distill sheets
	- This is mostly working however, it seems the genome_summary.xlsx output only has stuff from the ag sheet and not from the carbon sheet 
		- Need to fix this.
	- Alsom COMBINED_DISTILL() is happening twice
	- it seems the use of --distill_topic "default" is not capturing all of the default sheets
		- looks like in this case, carbon got put in twice

------

## Wednesday January 31
- Reworked distill_custom to only take in one file. 
	- This will require the user to cat all custom sheets together - not a bit ask, I feel.
- Seems grabbing ids from the combined_annotaions file is not working fully - testing shows dbcan are being pulled but not kofam
	- need to rename taget_id to kofam_id in the kofam_formatter
	- Then re-work the distill_summary to ensure it pulls all of the "\_id" column values
Tested command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-bin-qual-test --threads 15 --call --rename --annotate --use_dbcan --use_kofam --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/distill_ecosystems/test-2-topic.tsv -resume
```
- took quite a bit of tweaking but worked!
- distillate final output, `distillate.xlsx` consists of
	- individual sheets for each unique topic_ecosystem provided via the separate distillate sheets
			- by `--distill_topic <default|misc|energy|transport|carbon|nitrogen>` and/or `--distill_ecosystem <eng_sys|ag>` and/or `--distill_custom <path/to/TSV>`
	- individual sheets for tRNAs and rRNAs identified
		- **does not yet include flexibility of NOT providing these - need to rework channel IO and distill scripts**
	- DOES NOT include genome_stats. 
		- will be provided in part by the trna and rrna collected files, the combined_annotations file and a newly created combined_rrna file
		- DOES now - it works

------


## Thursday February 1
- Got distill_final working!
	- includes tRNA and rRNA sheets and provided distill_sheets
	- includes genome_stats with taxonomy, bin quality and rRNA
	- NEED to get in tRNA counts

------
## Friday February 2
- Distill.xlsx is complete! (for now)
	- includes tRNA and rRNA sheets and provided distill_sheets
	- includes genome_stats with taxonomy, bin quality and rRNA and tRNA counts
	- standalone rRNA and tRNA sheets

```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-bin-qual-test --threads 15 --call --rename --annotate --use_dbcan --use_kofam --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/test_data/distill_ecosystems/test-2-topic.tsv -resume
```



  - Next, going to focus on ensuring the apps can be run independently from one another. 
	- Annotate w/o call needed inputs:
		- protein fasta file (hmm searches) (.faa)
			- needs to have sample-prefixed headers (what we use --rename to rename fasta headers when calling genes)
		- nucleotide fasta file (mmseqs searches) (.fna)
			- needs to have sample-prefixed headers (what we use --rename to rename fasta headers when calling genes)
		- tRNAs (optional but cannot be run on called genes)
		- rRNAs (optional but cannot be run on called genes)
		- 

```
nextflow run DRAM2.nf --input_genes ../test_data/called_genes/ --input_proteins ../test_data/called_genes/ --outdir DRAM2-annot-only --threads 15 --annotate --use_dbcan --use_kofam --distill_topic 'carbon transport' --distill_ecosystem ag --taxa /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/HMP_DRAM1_test_data/gtdbtk.bac120.summary.tsv --rrnas /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/called_genes/rRNA/ --trnas /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/called_genes/tRNA/ -resume -with-timeline
```


Results:
- Worked! 
- However, the tRNA and rRNA column names need to be updated to new distillate column names
- Also, tRNA and rRNA sample naming comes from the file name from barrnap or tRNAscan-SE 
	- would have to require this
- Similarly, the input gtdb and checkm files have to have sample names that match
	- And I hardcoded swapping "." for "-" in the code.. would need to set requirements for this

Now testing on riviera:

```ruby
nextflow run DRAM2.nf --input_genes /test_data/called_genes/ --input_proteins /test_data/called_genes/ --outdir DRAM2-annot-only-Riviera --threads  --annotate --use_dbcan --use_kofam --distill_topic default --distill_ecosystem 'ag eng_sys' --taxa test_data/HMP_DRAM1_test_data/gtdbtk.bac120.summary.tsv --bin_quality test_data/HMP_DRAM1_test_data/checkm1-HMP-76-quality.tsv --rrnas test_data/called_genes/rRNA/ --trnas test_data/called_genes/tRNA/ -resume -with-timeline -with-report -with-dag
```
- Beefed up the masforks a ton! will have to check out the dag...

Results: worked! 

`/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/Testing-Results/DRAM2-annot-only-Riviera`

```ruby
rwoyda@l001:~/DRAM2-Nextflow/DRAM2-NF$ nextflow run DRAM2.nf --input_genes test_data/called_genes/ --input_proteins test_data/called_genes/ --outdir DRAM2-annot-only-Riviera --threads  --annotate --use_dbcan --use_kofam --distill_topic default --distill_ecosystem 'ag eng_sys' --taxa test_data/HMP_DRAM1_test_data/gtdbtk.bac120.summary.tsv --bin_quality test_data/HMP_DRAM1_test_data/checkm1-HMP-76-quality.tsv --rrnas test_data/called_genes/rRNA/ --trnas test_data/called_genes/tRNA/ -resume -with-timeline -with-report -with-dag -with-trace
Nextflow 23.10.1 is available - Please consider updating your version to it
N E X T F L O W  ~  version 23.10.0
Launching `DRAM2.nf` [silly_descartes] DSL2 - revision: b87d032505
WARN: Access to undefined parameter `use_pfam` -- Initialise it to a default value eg. `params.use_pfam = some_value`
WARN: Access to undefined parameter `use_merops` -- Initialise it to a default value eg. `params.use_merops = some_value`
WARN: Access to undefined parameter `distill` -- Initialise it to a default value eg. `params.distill = some_value`

DRAM2 Nextflow
===================================
fastas       : 0
outdir       : DRAM2-annot-only-Riviera
threads      : true
call genes   : false
annotate     : true
distill      : false

executor >  slurm (465)
[5e/1e0221] process > HMM_SEARCH_KOFAM (bin-18)    [100%] 76 of 76 ✔
[f8/f10450] process > PARSE_HMM_KOFAM (bin-18)     [100%] 76 of 76 ✔
[21/d635a4] process > KOFAM_HMM_FORMATTER (bin-18) [100%] 76 of 76 ✔
[8f/8d5f14] process > HMM_SEARCH_DBCAN (bin-26)    [100%] 76 of 76 ✔
[26/89e9fe] process > PARSE_HMM_DBCAN (bin-26)     [100%] 76 of 76 ✔
[78/414b1b] process > DBCAN_HMM_FORMATTER (bin-26) [100%] 76 of 76 ✔
[eb/f06ee5] process > COMBINE_ANNOTATIONS          [100%] 1 of 1 ✔
[fb/98ecdf] process > COUNT_ANNOTATIONS            [100%] 1 of 1 ✔
[eb/6ac585] process > ADD_BIN_QUALITY              [100%] 1 of 1 ✔
[d8/bf0c01] process > ADD_TAXA                     [100%] 1 of 1 ✔
[17/c4e86a] process > COMBINE_DISTILL (1)          [100%] 1 of 1 ✔
[ab/10be2b] process > DISTILL_SUMMARY (1)          [100%] 1 of 1 ✔
[18/a90ff6] process > TRNA_COLLECT                 [100%] 1 of 1 ✔
[a7/033a42] process > RRNA_COLLECT                 [100%] 1 of 1 ✔
[ab/c8fd93] process > DISTILL_FINAL (1)            [100%] 1 of 1 ✔
Completed at: 02-Feb-2024 18:52:46
Duration    : 16m 42s
CPU hours   : 23.8
Succeeded   : 465

```


Next, testing the 76 HMP bins on Riviera:
```groovy
nextflow run DRAM2.nf --input_fasta test_data/HMP_DRAM1_test_data/ --outdir DRAM2-annot-only-CAD-HMP-76-Riviera --rename --call Riviera --threads 20  --annotate --use_dbcan --use_kofam --distill_topic default --distill_ecosystem 'ag eng_sys' --taxa test_data/HMP_DRAM1_test_data/gtdbtk.bac120.summary.tsv --bin_quality test_data/HMP_DRAM1_test_data/checkm1-HMP-76-quality.tsv -resume -with-timeline -with-report -with-dag -with-trace
```

Worked!
Also re-vamped initial output:

```ruby
DRAM2 Nextflow
===================================
fastas       : test_data/HMP_DRAM1_test_data/
outdir       : DRAM2-annot-only-CAD-HMP-76-Riviera
threads      : 20
rename       : true
call genes   : true
annotate     : true
databases    : 
distill      : false
  topic      : default (carbon, energy, misc, nitrogen, transport)
  ecosystem  : ag eng_sys 
  custom     : none


executor >  slurm (5)
[5c/6638c7] process > RENAME_FASTA (bin-19)        [100%] 76 of 76, cached: 76 ✔
[22/0430ef] process > CALL_GENES (bin-26)          [100%] 76 of 76, cached: 76 ✔
[62/24c9ea] process > TRNA_SCAN (bin-19)           [100%] 76 of 76, cached: 76 ✔
[5c/50b4cc] process > TRNA_COLLECT                 [100%] 1 of 1, cached: 1 ✔
[47/1d1ca7] process > RRNA_SCAN (bin-29)           [100%] 76 of 76, cached: 76 ✔
[2a/b0e8c5] process > RRNA_COLLECT                 [100%] 1 of 1, cached: 1 ✔
[5f/8fcaaa] process > HMM_SEARCH_KOFAM (bin-66)    [100%] 76 of 76, cached: 76 ✔
[64/65fe4e] process > PARSE_HMM_KOFAM (bin-55)     [100%] 76 of 76, cached: 76 ✔
[ef/22e9cc] process > KOFAM_HMM_FORMATTER (bin-3)  [100%] 76 of 76, cached: 76 ✔
[e9/8953fc] process > HMM_SEARCH_DBCAN (bin-66)    [100%] 76 of 76, cached: 76 ✔
[31/b1468b] process > PARSE_HMM_DBCAN (bin-55)     [100%] 76 of 76, cached: 76 ✔
[b4/7d75c0] process > DBCAN_HMM_FORMATTER (bin-74) [100%] 76 of 76, cached: 76 ✔
[bf/1a556a] process > COMBINE_ANNOTATIONS          [100%] 1 of 1, cached: 1 ✔
[44/a99a5c] process > COUNT_ANNOTATIONS            [100%] 1 of 1 ✔
[2f/1c89b0] process > ADD_BIN_QUALITY              [100%] 1 of 1 ✔
[4d/1971c3] process > ADD_TAXA                     [100%] 1 of 1 ✔
[12/74207e] process > COMBINE_DISTILL (1)          [100%] 1 of 1, cached: 1 ✔
[93/0a8bf8] process > DISTILL_SUMMARY (1)          [100%] 1 of 1 ✔
[87/d01c96] process > DISTILL_FINAL (1)            [100%] 1 of 1 ✔

```


`/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/Testing-Results/DRAM2-annot-only-CAD-HMP-76-Riviera`

------

## Saturday February 3
- Annotate is now independent of call and can feed into distill
- Now, need to make distill independent
	- Distill relies on:
		- ch_combined_annotations - tsv of combined annotations (may or may not have taxonomy and bin quality metrics)
		- ch_bin_quality - user can provide with `--bin_quality` but is not required
			- user could also already have completeness and contamination already in the annots
			- Need to add in check for presence of completeness and contamination columns in combined_annots - if not present do not add to genome_stats
		- ch_taxa - tsv of GTDB taxonomy provided with `--taxa` but is not required
			- user could also already have taxonomy already in the annots
			- Need to add in check for presence of taxonomy column in combined_annots - if not present do not add to genome_stats
		- ch_collected_tRNAs - not required but individual tRNAs files from tRNAscan-SE can be provided with `--trnas`
			- Need to add in check for "null" contents in distill_final - if null then do not generate tRNA columns
		- ch_collected_rRNAs - not required but individual rRNAs files from barrnap can be provided with `--rrnas`
			- Need to add in check for "null" contents in distill_final - if null then do not generate rRNA columns

Testing command:
```groovy
# Using combined_annotations.tsv which has taxa and bin quality
nextflow run DRAM2.nf --annotations ../test_data/annotations/combined_annotations.tsv  --outdir DRAM2-distill-only-Riviera --threads 15 --distill_topic default --distill_ecosystem 'ag eng_sys' --rrnas ../test_data/called_genes/rRNA/ --trnas ../test_data/called_genes/tRNA/ -resume -with-timeline -with-report -with-dag
```
- Works!

```groovy
# Using combined_annotations_noBin-taxa.tsv which has neither taxa nor bin quality
# Providing these to be added
nextflow run DRAM2.nf --annotations ../test_data/annotations/combined_annotations_noBin-taxa.tsv  --outdir DRAM2-distill-only-Riviera --threads 15 --distill_topic default --distill_ecosystem 'ag eng_sys' --rrnas ../test-data/called_genes/rRNA/ --trnas ../test_data/called_genes/tRNA/ --taxa ../test-data/HMP_DRAM1_test_data/gtdbtk.bac120.summary.tsv --bin_quality ../test-data/HMP_DRAM1_test_data/checkm2-HMP-76-quality.tsv -resume -with-timeline -with-report -with-dag
```
- Works!

```groovy
# Using combined_annotations_noBin-taxa.tsv which has neither taxa nor bin quality
# Not providing these, dont add
nextflow run DRAM2.nf --annotations ../test_data/annotations/combined_annotations_noBin-taxa.tsv  --outdir DRAM2-distill-only-Riviera --threads 15 --distill_topic default --distill_ecosystem 'ag eng_sys' --rrnas ../test-data/called_genes/rRNA/ --trnas ../test_data/called_genes/tRNA/ -resume -with-timeline -with-report -with-dag

```
- Works!
  
  
  - Starting to work on MERGE_ANNOTATIONS()
	  - wrote process but need to create some fake annotations to test

------

## Wednesday February 7
- Working on merge_annotations()

```groovy
nextflow run DRAM2.nf --annotations ../test_data/annotations/combined_annotations_noBin-taxa.tsv --add_annotations ../test_data/annotations/combined_annotations_FAKE.csv --outdir DRAM2-merge-annots --distill_topic default --distill_ecosystem 'ag eng_sys' --rrnas ../test-data/called_genes/rRNA/ --trnas ../test_data/called_genes/tRNA/ -resume
```

----- 

## Friday February 9
- Working on merge_annotations() still
	- Looks like it is working!
	- Also added in the ability to provide EC numbers in a distill sheet.
		- wondering if this is desired for genbank... or should i just stop? :D 
However:
- **The current distill summary does not correctly output the additional columns. We need at add this back in eventually**
	-  as of now, additional columns in distill sheets gets a unique name in the output
	- also, all of these columns are in all final distill output sheets

-----

## Monday February 12
- Taking a step back from the main apps and going to implement some more databases
	- Starting with vogdb (hmmsearch)
		- created if statement and logic for vog database searching
		- testing hmmsearch and parse_hmm
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-VOGdb-test1-smaller --threads 2 --call --rename --annotate --use_vog --use_dbcan --distill_topic default --distill_ecosystem 'eng_sys ag' -resume 
```

Results: worked! Woot!


- [ ] Working on CAMPER next - processes written - need to verify test worked:

Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-VOGdb-test1-smaller --threads 2 --call --rename --annotate --use_vog --use_dbcan --use_camper --distill_topic default --distill_ecosystem 'eng_sys ag' -resume
```

Results:
- Worked!
- However, the distill camper sheet has many additional columns and highlights the need for separating out the additional columns in the distillate.xlsx


-----


## Tuesday February 13th
- Going to start working on the mmseqs2 blast style search-databases
- Starting with MEROPS
	- Wondering if I should check if a mmseqs2 searhc is going to occur and then index all the inputs
		- Lets just index separately

Testing MMSEQS_INDEX():
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-MEROPS-test --threads 2 --call --rename --annotate --use_merops
```
- Works!
Bigger test:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-MEROPS-test --threads 2 --call --rename --annotate --use_merops --use_vog --use_dbcan --use_camper  --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/distill_camper_Feb122024.tsv -resume
```

Results:
WORKED! Woot!

Next, work on the renaming mmseqs databases (recalling that pfam is a bit different)
- [x] refseq_viral
- [ ] kegg
- [ ] uniref
- [ ] pfam

------

## Wednesday February 14

Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-MEROPS-test --threads 2 --call --rename --annotate --use_viral
```

- Works
Bigger test:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-MEROPS-test --threads 2 --call --rename --annotate --use_merops --use_viral --use_dbcan --use_camper  --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/distill_camper_Feb122024.tsv -resume
```
- Works!

- Also, was able to get away from combining the collected_formatted_hits in the way I was - each db as its own .mix()
	- Now, it will mix these when the annotation has occured and then does a final mix to create one output channel

Now, get these database implemented:
- [x] camper mmseqs
	- [ ] test
- [ ] canthyd
	- [x] mmseqs
		- [ ] test
	- [ ] hmm
- [ ] fegenie hmm
- [x] methyl mmseqs
	- [ ] test
- [ ] sulfur hmm

Some of the MMseqs databases above also come with a list which needs to be incorporated into the formatted output
- Really just description and rank info

Before I do this I will test that the above mmseqs work without their associated list(s):
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-MEROPS-test --threads 2 --call --rename --annotate --use_cant_hyd --use_methyl --use_camper
```
- works!


## Thursday February 15

Also, I added in functionality to the combine_annoations to keep query_ids which have unique start and stop - beofre this was overwriting start and stop and thus only keeping one start and stop even if there were multiple annotations for a query_id


Test on all databases so far:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-all-db-test --threads 2 --call --rename --annotate --use_cant_hyd --use_camper --use_dbcan --use_vog --use_viral --use_methyl --use_merops --distill_topic default --distill_ecosystem 'ag eng_sys' --distill_custom assets/forms/distill_sheets/distill_camper_Feb122024.tsv -resume
```

- OMG it worked!
  
  Next, going to migrate over to Riviera and do a big test on all of the 76 HMP samples with the command above.
  
  Command:
  ```groovy
nextflow run DRAM2.nf --input_fasta test_data/HMP_DRAM1_test_data/ --outdir DRAM2-all-db-CAD-HMP-76-Riviera --threads 20 --call --rename --annotate --use_cant_hyd --use_camper --use_dbcan --use_vog --use_viral --use_methyl --use_merops --distill_topic default --distill_ecosystem 'ag eng_sys' --distill_custom assets/forms/distill_sheets/distill_camper_Feb122024.tsv -resume
```

Results: DONE - Worked!
- Had to change all relative paths in the nextflow.config to absolute paths
	- may need to deal with this later...
- Annotations
	- in general is correct
	- the order of the columns needs to be dealt with but it is not too far off
	- distill needs to dedupe the gene_id column
		- This is occurring now because we introduced dupes when we updated the start_position and stop_position 


Next, going to keep adding mmseqs databases. 
- [x] KEGG 
	- [ ] Needs work
- [x] pfam
	- [ ] Needs mmseqs db rebuilding
- [ ] uniref
	- [ ] not tested yet

KEGG test:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-kegg-pfam-uniref-db-test --threads 10 --call --rename --annotate --use_kegg -resume
```
- Worked!
- Need to dedupe distillate.xlsx
- The taxa and bin quality is not present on the genome_stats - look into this
	- The taxa and bin quality info is in the annotations
- However, i just now descovered what the SQL database is for
	- This holds all of the info to match descriptions to gene IDs
	- 
		- Ex: the KEGG output gene_ids are not KOs they are kegg gene IDs and to get KOs we need to 
			- Option 1: Use the SQL database 
				- Not a huge ask, and would be nice to have all the descriptions in one place
				- But, this database is 36G and would need to be provided as well (or built from the user's end but this just adds more room for user-errors)
			- Option 2: Use spreadsheets like I have been doing
				- This is not bad but I want the best efficiency and I am not sure this is better than spreadsheet parsing - going to ask Alex
				- This means for KEGG, for example, that I would need to parse the .pep files and produce a 4-column spreadsheet: gene_id|KO|description|EC

------
## Friday February 16
- Starting with testing pfam and then uniref


Pfam:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-kegg-pfam-uniref-db-test --threads 10 --call --rename --annotate --use_pfam --distill_topic default --distill_ecosystem 'ag eng_sys' -resume
```

Results:
- Failed
- This is due to pfma being a mmseqs profile and not a mmsdb file
	- Need to determine why pfam was done this way and if I can change it to a "normal" mmsdb file
	- 

------

## Monday February 26

- Need to do 2 major things this week and 1 minor
	- Get the SQL database working
		- Just using the existing SQL database from: `give/path`
	- Get distill and working flawlessly
	- Get pfam mmseqs db reformatted and working

### SQL database
- using the descriptions_db.sqlite from DRAM 1.4: `/home/Database/DRAM/sep_12_22_dram1.4.0_full_db/description_db.sqlite`
- wrote nextflow process and supporting python script to match gene_id values and add the description

Testing on uniref:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-uniref-db-test --threads 20 --call --rename --annotate --use_uniref --distill_topic default --distill_ecosystem 'ag eng_sys' -resume --slurm_node zenith
```

Results:
Worked!

Next, add this in for all annotators which have a database table:
```
dbcan_description      pfam_description       vogdb_description    
kegg_description       uniref_description   
peptidase_description  viral_description
```
I renamed peptidase in the sql database:
```
dbcan_description   merops_description  uniref_description  vogdb_description 
kegg_description    pfam_description    viral_description
```

Adding to the other mmseqs searches which have sql description tables:
- VIRAL
- MEROPS
- KEGG
  
  
Testing:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/subsampled/test_2/smaller/ --outdir DRAM2-uniref-db-test --threads 10 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --distill_topic default --distill_ecosystem 'ag eng_sys' -resume --slurm_node zenith
```
Results:

- KEGG
	- Needs additional processing to pull out the KO - it is in parentheses in the description
		- Either pass a flag to the existing process or create a different process which can do this
		- DONE!
- MEROPS
	- Worked!
	- Does not have EC
- VIRAL
	- Worked!
	- Does not have EC

Thus, some databases have extra info we need. 
- Adding in checks for kegg and (eventually DBCAN) 


Bigger test with taxa, bin_quality and RNAs on all databases (pfam still not working)
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 10 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog -resume --slurm_node zenith
```
Results:
Running...


- Added in SLURM job names for each process. 
	- This is only to make the squeue output more informative - now they all start with "DRAM2-..."


Next, need to add this functionality for the HMM searches:
Databases:
 - [x] KOFAM (NOT NEEDED)
 - [x] DBCAN
	 - Testing - ...
 - [x] CAMPER (NOT NEEDED!)
 - [x] VOGdb


------

## Tuesday February 27
- Working on adding a ${db_name} == dbcan check to the add_sql_descriptions.py script 
	- Removed adding info from the dbcan hmm formatter
- Got dbcan to work! 
- next, vogdb
	- Worked!
- Next, cant_hyd hmm 
	- 
Testing canthyd hmm:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 10 --call --rename --annotate --use_canthyd -resume --slurm_node tardigrade
```
Results: works!
- Had to do some tinkering. Also got the rank working - Need to do this for other datbases with a rank


Next, add this same rank functionality to the CAMPER hmm search results
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 10 --call --rename --annotate --use_camper -resume --slurm_node tardigrade
```
Results:
- Worked!

Next database, sulfur:
Testing sulfur hmm:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 10 --call --rename --annotate --use_canthyd -resume --slurm_node tardigrade
```
Results:
- Worked!

Next, pfam HMM (i forgot about this one because the mmseqs is being tricky)
**This does not need to be done. I was wrong"
- The pfam hmm is only a descriptions file and I think i will forego uising it and rely on the info in the sql database


Next, fegenie HMM:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-fegenie --threads 10 --call --rename --annotate --use_fegenie -resume --slurm_node tardigrade
```
Results: Worked!


**As of here, all databases are implemented and working except for PFAM mmseqs**

Testing all dbs again:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 10 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur -resume --slurm_node tardigrade
```
Worked!

Will run bigger test on Riviera soon.


Messing with max forks:
    max_forks_single_cpu  = 2
    max_forks_user_cpu = 1
    max_forks_hmm = 2
- Where max_forks_single_cpus is used for processes which only utilize 1 cpu - thus, more of them can run
- And max_forks_user_cpus is used for processes which utilize the number of cpus the user specified - thus, more of them can run

Command to test:
with 
    max_forks_single_cpu  = 2
    max_forks_user_cpu = 1
    max_forks_hmm = 2
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_10/ --outdir DRAM2-HMP-10-test-max-forks --threads 10 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --slurm_node tardigrade -with-report -with-timeline
```

- This was taking too long with Uniref and KEGG
	- Going to take these out of the test
		- Will test a BIG run on Riviera once it is back up

-------

## Wednesday February 28
- Going to mess with distill today


Need to see where its at, will run:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 10 --call --rename --annotate --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default -resume --slurm_node tardigrade
```
Results:
- distill summary failed


------
## Friday March 1
- Got distill working!

Command tested:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/test_2/ --outdir DRAM2-HMP-2-test-all-dbs --threads 5 --call --rename --annotate --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_kofam --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```


However, the EC partial matching is not working
- I did get it to work but the sheer number of possible matches with the EC hierarchy can be too much
	- That is, if the user says get me all 3.- distill numbers we can easily output over the max allowed line for a xlsx document.
- For now, exact matches are done for both gene_ids and ec_ids


Next, going to run some larger-scale test but on subsampled inputs.
- We will use the HMP 76 but subsample them
- Use all databases - will take a bit. 

Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/76_subsampled/ --outdir DRAM2-HMP-76-subsampled-12-dbs-CAD-03012024 --threads 5 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith -with-report -with-trace -with-timeline
```

Results:
- It worked! 
- took 6hrs
- Distill all looks great except for
	- unwanted extra "gene_id" column in the EC_test distill sheet
	- need to re-order the genome stats output columns
- In general, I am very pleased!\

------

## Monday March 4
- Updating help menu

------


## Wednesday March 6 and Thursday March 7
- Fixing gene start and stop positions
	- These positions, which are in the output `raw-annotations.tsv`, are currently the annotation hit start and stop on the contig
	- Needs to just be the called gene start and stop position
	- Adding in a tsv generator for a gene_locs.tsv
		- This is passed into both the mmseqs and hmm formatter scripts and sets the gene start and stop positions
- Problem: ran into the issue I did with COMET when passing in two sample tuples into a process
	- These need to be combined first and then sent into the process to ensure they have the same sample value

Here is code from COMET:
```groovy
            // Combine queue channels into one
            checkm_input = MAP_BIN.out.metabat_bins.join(MAP_BIN.out.metabat_bin_depths).map{ item -> 
                    [item[0], item[1], item[2], item[3], item[6]]}
```


- Also, need to add in functionality to produce the gene_locs tsv file if the user provides input genes, which can be a .gff or .fna,
	- Then, this can be combined into the correct channel for mmseqs search and hmm search formatters

Added in functionlisty for combining the channels but need to test


Testing command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-methyl-gene-locs-test-03062024 --threads 5 --call --rename --annotate --use_methyl -resume --slurm_node zenith
```
- Works!

I altered this for all mmseqs searchs but needs to be done for each hmm formatter:
- [x] kofam
- [x] dbcan
- [x] camper
- [x] fegenie
- [x] canthyd
- [x] sulfur 
- [ ] vogdb



- Modifying add_annotations
- Then, working on merge_annotations


add_annotations:
Test command(s):

Generate initial annotations which we will use to add: (kofam and dbcan)
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-kofam-dbcan-03062024 --threads 5 --call --rename --annotate --use_kofam --use_dbcan --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```
- DONE

Starting at Call: (add in kofam and dbcan to all other dbs (excluding pfam))
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-add-annots-CAD-03062024 --threads 5 --call --rename --annotate --add_annotations DRAM2-HMP-2-subsampled-kofam-dbcan-03062024/RAW/raw-annotations.tsv --use_merops --use_viral --use_camper --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```
- Worked!

Starting at Annotate: (add in kofam and dbcanto all other dbs (excluding pfam))
```groovy
nextflow run DRAM2.nf --input_genes DRAM2-HMP-2-subsampled-add-annots-CAD-03062024/Prodigal_v2.6.3/ --outdir DRAM2-HMP-76-subsampled-add-annots-AD-03072024 --threads 5 --annotate --add_annotations DRAM2-HMP-2-subsampled-kofam-dbcan-03062024/RAW/raw-annotations.tsv --use_merops --use_viral --use_camper --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-1-multiple-gene-id.tsv --slurm_node zenith -resume --taxa ../test_data/HMP_DRAM1_test_data/gtdbtk.bac120.summary.tsv --bin_quality ../test_data/HMP_DRAM1_test_data/checkm2-HMP-76-quality.tsv --trnas DRAM2-HMP-2-subsampled-kofam-dbcan-03062024/tRNA_rRNA/collected_trnas.tsv --rrnas DRAM2-HMP-2-subsampled-kofam-dbcan-03062024/tRNA_rRNA/collected_rrnas.tsv
```
- For starting at annotate and not call, we run into the issue of not having the gene locations
- If we require the user to give a .fna FASTA file, we can technically take the locations out of here and generate a table
	- However, this requires the user to make sure they have a FASTA with the gene locations 
- [x] Added in GENE_LOCS() process which will run when `--call` is not provided
	- [x] assumes input `--input_genes *.faa` files with gene start and stop locs in positions 2 and 3 in the header


For starting at distill, need the base annotations:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-10dbs-CAD-03062024 --threads 5 --call --rename --annotate --use_merops --use_viral --use_camper --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```

Starting at Distill:
```groovy
nextflow run DRAM2.nf --outdir DRAM2-HMP-76-subsampled-add-annots-D-03062024 --threads 5 --annotations DRAM2-HMP-2-subsampled-10dbs-CAD-03062024/raw-annotations.tsv --add_annotations DRAM2-HMP-2-subsampled-kofam-dbcan-03062024/RAW/raw-annotations.tsv --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith 
```

------


## Monday March 11
- Still need to finish --add_annotations testing from March 6th and 7th
	- Last steps
- [x] Rebuild main singularity image with QUAST
	- [x] Build QUAST process
		- [x] combine results into single CSV and pass to distill
Testing command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-QUAST-test-03112024 --threads 5 --call --rename --annotate --use_camper --use_dbcan --distill_topic 'carbon transport camper' --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```
- [x] Test the main container with --rename

- [x] alter mmseqs search to have pfam and kegg checks
	- [x] if KEGG then do reverse best hits and additional filtering
		- [ ] Need to test
	- [ ] if pfam then do the profile search
		- [x] This is being tricky
			- [x] The alignment keeps failing - updated mmseqs2 to the newest version - works
				- [x] I was able to run a super small test outside of DRAM2 which worked - not sure why is failing inside DRAM2
				- [ ] Doing:
					- [x] re-downloading pfam fasta and will rebuild the mmseqs2 profile - Will take all day - test tomorrow if new image does not work
					- [x] re-building singularity image with newest version of mmseqs2 - there is one more version ahead of what we were using
						- [x] Testing.. - OMG it worked! 
						- [ ] Double check output


- [ ] Re-organize annotations output
	- [ ] put kegg first if present
	- [ ] remove strandeness column
	- [ ] group database columns together

Testing command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsam-kegg-annot-test-03112024 --threads 10 --call --rename --annotate --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd -resume --slurm_node zenith
```
- had to cancel will resume tomorrow


------


## Tuesday March 12
- Focusing on testing KEGG RBH and the re-ordered `raw-annotations.tsv`
- Next, will go back to March 6th and 7th and finish testing `--add_annotations`
- Then, work on `--merge_annotations`

- KEGG RBH and `raw-annotations.tsv` re-ordering
Testing command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsam-kegg-annot-test-03112024 --threads 10 --call --rename --annotate --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd -resume --slurm_node zenith
```
- KEGG reverse searching takes forever!
	- I commented this out for now. will test later on Riviera
		- Going to re-run and see if the raw annotations reordering works
- DONE - works!


Add_annotations testing: (taken from March 6th and 7th)
REDO: Generate initial annotations which we will use to add: (kofam and dbcan)
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-kofam-dbcan-03062024 --threads 5 --call --rename --annotate --use_kofam --use_dbcan --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```

For starting at distill, need the base annotations:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsampled-8dbs-CAD-03122024 --threads 5 --call --rename --annotate --use_merops --use_viral --use_camper --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith
```
- Worked!

Starting at Distill:
```groovy
nextflow run DRAM2.nf --outdir DRAM2-HMP-2-subsampled-add-annots-D-03122024 --threads 5 --annotations DRAM2-HMP-2-subsampled-8dbs-CAD-03122024/RAW/raw-annotations.tsv --add_annotations DRAM2-HMP-2-subsampled-kofam-dbcan-03062024/RAW/raw-annotations.tsv --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv -resume --slurm_node zenith 
```
- Worked!

Finishing touches on formatting EC numbers:
- [x] dbcan
- [x] kegg
- [x] camper
	- [x] hmm
	- [x] mmseqs - NOT NEEDED

Testing command for dbcan, kegg and camper EC formatting:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsam-kegg-camper-dbcan-EC-test-03112024 --threads 20 --call --rename --annotate --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd -resume --slurm_node zenith
```
Results:
- Worked!
`DRAM2-HMP-2-subsam-kegg-camper-dbcan-EC-test-03112024`


Adding Rank:
For now we will use this scheme and update it later:
```
def assign_rank(row): # Dynamically check for database hits and assign ranks based on bitScore criteria databases = ["kegg", "uniref", "pfam", "dbcan", "merops", "vogdb"] db_scores = {db: row.get(f"{db}_bitScore", None) for db in databases} # Rank logic if db_scores["kegg"] is not None and db_scores["kegg"] > 350: return 'A' elif db_scores["uniref"] is not None and db_scores["uniref"] > 350: return 'B' elif (db_scores["kegg"] is not None and db_scores["kegg"] > 60) or (db_scores["uniref"] is not None and db_scores["uniref"] > 60): return 'C' elif any(db in row for db in ["pfam_id", "dbcan_id", "merops_id"]) and all(score <= 60 for score in db_scores.values() if score is not None): return 'D' elif db_scores["vogdb"] is not None or all(score is None or score < 60 for score in db_scores.values()): return 'E' return 'E' # Default to rank E if no other conditions are met
```


Testing:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsam-kegg-camper-dbcan-EC-test-03112024 --threads 20 --call --rename --annotate --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd -resume --slurm_node zenith
```
Results:
- Working! 

Generate test data for Maddie:
Command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-subsampled-13-dbs-CAD-03122024 --threads 20 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic 'default camper' --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv --slurm_node tardigrade -with-report -with-trace -with-timeline
```
Started on March 12th at 4:30 PM

Results:
- This is still running at 8:30AM on Wednesday
- Need to think about more efficient ways to do combine annotations - this is currently at 3 hours of processing


-------


## Wednesday March 13
- Working on a parallel version of combine_annotations() because the HMP-76 dataset with 13 databases was up to 3 hours.

Testing parallel version of combine_annotations:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsam-kegg-camper-dbcan-EC-test-03112024 --threads 20 --call --rename --annotate --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd -resume --slurm_node zenith
```
- Worked!
  
Going to resume the HMP 76:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-subsampled-13-dbs-CAD-03122024 --threads 20 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic 'default camper' --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/test-3-EC.tsv --slurm_node tardigrade -with-report -with-trace -with-timeline -resume
```
Results:
- combine annotations ran in like 15 min (?)!


- I also need to run the HMP 76 dataset through DRAM1 and give to Maddie:
Script:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=reedrich@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --nodelist=tardigrade


# activate the DRAM2 environment
source /opt/Miniconda2/miniconda2/bin/activate DRAM2BETA

dram2 -d /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/rory_DRAM2_HMP-76_analysis -t 30 call /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/HMP_DRAM1_test_data/*.fa

dram2 -d /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/rory_DRAM2_HMP-76_analysis -t 30 annotate --use_dbset adjectives_kegg_set --use_db camper

dram2 -d /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/rory_DRAM2_HMP-76_analysis -t 30 distill


```
Started at 9:03AM
- Want to compare this to the total time for DRAM2-NF (however, i have no idea how to compare this b/c of threads and max_forks..)




What else to work on today:
- [ ] Partial EC searching - pick up on this tomorrow
	- [ ] correct counting in distill output cols
- [ ] Correct distill output and counting when multiple genes given in gene_id distill input sheet column (example: CAMPER)
	- [ ] correct counting in distill output cols
- [ ] Generate GFF and GBK from annotations
	- [x] GFF - initial try is working!
		- [ ] Should give this to Henry at KBase team to look at
			- [ ] Still need to implement a user-based list of databases to include.
	- Need biopython for gbk generation
		- [x] Re-building current version of the container
			- [ ] GBK initial try is working
				- [ ] pass to Kayla for inspection
- [ ] MIMAG and MISAG column in Genome_Stats page
	- [ ] implement logic in distill_xlsx.py




Going to resume the HMP 76 run again overnight and give results to Maddie:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-subsampled-13-dbs-CAD-03122024 --threads 20 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -with-report -with-trace -with-timeline -resume
```
Results:
- Worked!
- Going over results with Kayla at 10 AM March 14th
  
  
--------

## Thursday March 14
### Recap:
- GFF and GBK are generated but will need touching up after sending them to Kayla, B, and Henry
- Partial EC matching and counting needs work
	- I need to step back and rethink how this needs to work and what the output distillate counts should be reflective of
- Need MIMAG output column in distillate.xlsx


### What next:
- [x] Independent Merge annotations
- [x] Independent fasta header rename option
- [x] Partial EC matching and counting
- [ ] Multi EC counting
- [ ] MIMAG output distill column

### Progress
- [x] Merge annotations
	- [x] first generate a bunch if different test annotation files

Test command:
```groovy
nextflow run DRAM2.nf --outdir DRAM2-merge-annots-test-3-03142024 --threads 20 --merge_annotations ../test_data/merge_annots_test/--slurm_node tardigrade
```
- working!

- [x] Independent RENAME:
Testing command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-rename-test-03142024 --rename
```
- works!


### Meeting with Kayla 03142024
- use the 152 Gary's MAGs and jareds 10 mags as test data for Maddie test data
	- `/home/projects-wrighton-2/Permafrost_18O/results/metabat_v2.15/all_metabat_bins`
		- I have put all of these bins here: `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test`
	- can still give her the 76 HMP output
- [ ] update query_id sorting
	- [ ] need to add prefix 0s so that they sort correctly
- test tRNA and rRNA identification on 250 bp and 2500 bp




Small test on Gary + Jared MAGs:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/test_2/ --outdir DRAM2-Gary-Jared-2-test-CAD-03142024 --threads 20 --call --rename --annotate --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume
```
- [x] Need to change rRNA scan to ensure even if we find nothing that we output a tsv with "NULL" in it. 
	- [x] Also had to generalize add tax and bin quality processes as they had the gtdb and checmk1/2 hardcoded first col names. no it just assumes the sample names are in the first column
	- [ ] This is still broken and needs to be fixed
- [x] Also had to update the sql annotations database creation to account for the rank and geneo no. columns
- worked!

Plan to run the whole set over night:
- need the checkm data first.. and toxonomy.. may need to wait for these results to finish
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/ --outdir DRAM2-Gary-Jared-ALL-CAD-03142024 --threads 20 --call --rename --annotate --use_kegg --use_kofam --use_vog --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume
```  
Results:
Will have to start this once I have ran gtdb and checkm2

Alternative, run 76 HMP again overnight:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-subsampled-13-dbs-CAD-03142024 --threads 20 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -with-report -with-trace -with-timeline -resume
```
Results:
- This worked however, there is a bug for kegg annotations
	- In the sql add descriptions script we need to not produce a kegg_orthology column because in the future we are looking for gene ids in the _id-named columns
	- Thus, I have replaced the data in the kegg_id column and replaced it with the kegg_orthology values.
	- This is okay because the preexisting kegg_id values are actually contained within the kegg_description values if the user needs these.
- Re-running morning of March 15th with this change
Results:
- Worked!
- Got a ton more annotations in the distill!

--------

## Friday March 15
### Thursday recap:
- discovered some bugs when trying to run Gary and Jared's MAGs
- if no samples have rRNA identified, distill will crash because there are headers in these "empty" files
	- need to set rRNA file to "NULL" if only a header
- gene no. was causing the SQL database to fail
	- switched to gene_number
	- This count is wonky and needs to be corrected (see example below)
- got the partial EC matching to work
	- still needs a big test
- got merge_annotations to work


### What next:
- [x] Rerun the 76 HMP test from last night with kegg_id/kegg_orthology change
- [x] Need to re-do how the gene no. column is numbered.. there are some inconsistencies (see example below)
- [ ] MIMAG standards calculation and column in the distill genome stats page
- [x] Set rRNA and tRNA collected output files to "NULL" if nothing was identified
- [x] Generate Gary and Jared test dataset and run on Rory's DRAM2 (to get product) and run on DRAM2-NF to test it
	- [x] GTDB and checkm2 input file (will just be a single file with both GTDB and CHECKM2 information) - GTDB running now
- [ ] Optimize distill to make it faster


### Progress:
- [x] Need to fix gene_number counts in raw-annotations.tsv
Example: These are not in order w.r.t. gene start and stop locations

|                                            |                            |      |      |     |     |     |
| ------------------------------------------ | -------------------------- | ---- | ---- | --- | --- | --- |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_1  | blaze10-rt-cf-qt_MH_bin-10 | 2    | 403  |     | E   | 1   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_10 | blaze10-rt-cf-qt_MH_bin-10 | 7329 | 7841 |     | D   | 2   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_11 | blaze10-rt-cf-qt_MH_bin-10 | 8162 | 9319 |     | D   | 3   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_2  | blaze10-rt-cf-qt_MH_bin-10 | 403  | 1311 |     | D   | 4   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_3  | blaze10-rt-cf-qt_MH_bin-10 | 1311 | 2303 |     | E   | 5   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_4  | blaze10-rt-cf-qt_MH_bin-10 | 2554 | 3456 |     | D   | 6   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_5  | blaze10-rt-cf-qt_MH_bin-10 | 3565 | 4047 |     | D   | 7   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_6  | blaze10-rt-cf-qt_MH_bin-10 | 4312 | 5586 |     | D   | 8   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_7  | blaze10-rt-cf-qt_MH_bin-10 | 5689 | 6165 |     | D   | 9   |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_8  | blaze10-rt-cf-qt_MH_bin-10 | 6268 | 6843 |     | D   | 10  |
| blaze10-rt-cf-qt_MH_bin-10_k121_1103059_9  | blaze10-rt-cf-qt_MH_bin-10 | 6824 | 7291 |     | D   | 11  |
Testing command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/test_2/ --outdir DRAM2-Gary-Jared-2-test-CAD-03152024 --threads 20 --call --rename --annotate --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume
```
Results:
- works!


- [x] Also, fix when no rRNAs or tRNAs are found
	- [x] This means setting stuff = "NULL" and checking for it
	  - [x] tRNA
	  - [x] rRNA

Testing command for modifications to rRNA scripts:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/test_2/ --outdir DRAM2-Gary-Jared-2-test-CAD-03152024 --threads 20 --call --rename --annotate --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume
```
- This is a good test dataset b/c both of the inputs do not have any rRNAs
Results:
- Worked!

Testing command for modifications to tRNA scripts:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/test_2/ --outdir DRAM2-Gary-Jared-2-test-CAD-03152024 --threads 20 --call --rename --annotate --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume
```
- This is a not a good test set b/c no samples lack tRNAs, however we can test it to ensure it still works
Results:
- Worked!


Distill multiprocessing testing command:
This took about 1 hour previously to distill:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-subsampled-13-dbs-CAD-03142024 --threads 20 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -with-report -with-trace -with-timeline -resume
```
Results:
- Still took over an hour to distill
- Need to revisit the optimization implemented
- The rrna_sheet.tsv was not incorporated correctly because it was not generated correctly. See the distill working dir:
  `less work/9e/90f386128ecce954ed35ab7f7872a1/rrna_sheet.tsv`
  - [ ] Need to look into this
	  - [ ] I modified the code and save it. have not git pushed it. 
		  - [ ] have to wait until the Gary and Jared test is finished below

DRAM2 Big test: Gary and Jared's MAGs
data in: `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test`
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/ --outdir DRAM2-Gary-Jared-162-CAD-03152024 --threads 20 --call --rename --annotate --use_kegg --use_uniref --use_kofam --use_vog --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume
```  
- Started at 3:45PM on Friday
- Still running at 7AM Saturday
	- Finished! All worked execpt for:
		- rRNA needs to be fixed
		- I manually put in bin quality and taxa info b/c i forgot to remove the working dir version of the combined annotations
	- Going to re-resume after I test the rRNA fix on a smaller test data set


Generate Rory's DRAM2 Gary and Jared's MAGs output:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=reedrich@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --nodelist=tardigrade

# activate the DRAM2 environment
source /opt/Miniconda2/miniconda2/bin/activate DRAM2BETA

dram2 -d /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test/rory_DRAM2_Gary_Jared_162_analysis -t 20 call /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test/*.f*

dram2 -d /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test/rory_DRAM2_Gary_Jared_162_analysis -t 20 annotate --use_dbset adjectives_kegg_set --use_db camper

dram2 -d /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test/rory_DRAM2_Gary_Jared_162_analysis -t 20 distill
```
- This was a big fail.. not exactly sure why...
  
Re-running with DRAM1:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=400gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=reedrich@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low
#SBATCH --nodelist=tardigrade

# activate the DRAM environment
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4

DRAM.py annotate -i '/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/gary_jared_test/*.f*' -o DRAM1_Gary_Jared_162_analysis --use_camper --use_fegenie --use_sulfur --threads 20

DRAM.py distill -i DRAM1_Gary_Jared_162_analysis/annotations.tsv -o DRAM1_Gary_Jared_162_analysis --threads 20
```
- Started Monday at 8:45 AM
- Finished on Tuesday around 11:00AM



------

## Monday March 18
### Recap:
- Gary + Jared
	- DRAM2 (rory) failed hard
		- Re-running with DRAM1
	- DRAM2-NF worked
		- Need to fix rRNA issue and re-run
		- Also, remove the old combined annotations to fix bin quality and taxa incorporation


rRNA fix and test:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-subsam-kegg-camper-dbcan-EC-test-03112024 --threads 20 --call --rename --annotate --use_kegg --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd -resume --slurm_node tardigrade
```
- Finally got it to work

Now, re-run the Gary + Jared set:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/ --outdir DRAM2-Gary-Jared-162-CAD-03152024 --threads 20 --call --rename --annotate --use_kegg --use_uniref --use_kofam --use_vog --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume --taxa ../test_data/gary_jared_test/gary_jared_taxa_binqual_162.tsv --bin_quality ../test_data/gary_jared_test/gary_jared_taxa_binqual_162.tsv
```  
Results:
- Running...
	- rRNA worked!
	- add taxa and add bin quality did not work... because I forgot to add the flag!! (edited command above)
		- Restarted...
	- tax and bin quality worked!
	- Distill still running..
	- Dang - just realized the kegg output columns in the annotations are messed up. Going to modify `add_sql_descriptions.py`
	- Edited and restarted...
		- Worked!
		- Still waiting on distill....
- Done!
	- distill took forever.
		- [ ] look up how long 
	- [ ] Send location of data to Maddie
		- `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/DRAM2-Gary-Jared-162-CAD-03182024`
		- Also send DRAM1 output
			- [ ] Still running...

Run 76 HMP again overnight:
```groovy
nextflow run -bg DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-subsampled-13-dbs-CAD-03182024 --threads 8 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -with-report -with-trace -with-timeline
```
Results:
- Started at 8:09 PM Monday
- Finished
- 


--------


## Wednesday March 20
- [ ] Created todo list w.r.t. the march 26th deadline: [[March 20 2024 Daily Log Wednesday]]


### Recap
- DRAM2 is in a good place and the raw-annotations.tsv and distillate.xslx output is stable
- Still need minor tweaks to distillate.xlsx
	- add in the MIMAG categorization column to genome stats page
	- also need work on distill process optimization
- Still need to tweak GFF and Genbank output format
- GitHub needs a lot of attention



### Today
- [x] GFF and GBK formatting


GFF and GBK formatting test command:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-kofam-dbcan-GFF-GBK-03202024 --threads 20 --call --rename --annotate --use_kofam --use_dbcan --generate_gff --generate_gbk -resume --slurm_node zenith
```

- GFF is solid
	- [ ] pass to Henry and Kayla for review
	- [ ] bigger tests

- GBK is solid as well
	- [ ] pass to Henry and Kayla for review
	- [ ] bigger tests

- [ ] Need to put in the documentation what the inputs tables for taxa and bin quality have to be
	- [ ] recall, this is based on checkm and gtdb so we will require their format and column names
	- [ ] "user_genome", and "classification" for gtdb and "Bin Id", "Completeness" and "Contamination" for checkm1 and 2




## Thursday March 21
- [x] Trna null setting if all input trnas, to trna collect, are empty, null.
- [ ] reverse best hits for KEGG testing - just need to run a test to see how long this actually takes
- [ ] Distill optimization
- [x] Change custom distill to be directory based, not a single-file. 
	- [x] This will result in the topic and ecosystem coming into combine distill and the custom coming in as a list of paths
- [x] Create a default product process and get it set up for Maddie
- [ ] Start formulating some large tests
	- [ ] Create diverse test set: `/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/HMP_Gary_Jared_others`
		- [ ] local machine files: `/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/DRAM2-NF-test-datasets`
		- [ ] sample names with diverse characters
		- [ ] input MAGs and assemblies (metagenomic and single bacteria)
			- [ ] with diverse headers
	- [ ] Add annotations test
	- [ ] Merge annotations test
	- [ ] Rename test test
	- [ ] Call independent test
	- [ ] Annotate independent tests
	- [ ] Distill independent tests


Testing tRNA_collect update:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-kofam-dbcan-GFF-GBK-03212024 --threads 10 --call --rename --annotate --use_kofam --use_dbcan --generate_gff --generate_gbk --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/custom_distill/ -resume --slurm_node zenith
```
- works!

Testing distill custom directory:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-kofam-dbcan-GFF-GBK-03212024 --threads 10 --call --rename --annotate --use_kofam --use_dbcan --generate_gff --generate_gbk --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom assets/forms/distill_sheets/custom_distill/ -resume --slurm_node zenith
```
- works!

Testing distill custom empty directory:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-kofam-dbcan-GFF-GBK-03212024 --threads 10 --call --rename --annotate --use_kofam --use_dbcan --generate_gff --generate_gbk --distill_topic default --distill_ecosystem 'eng_sys ag' --distill_custom . -resume --slurm_node zenith
```
- works!
	- gives correct error: `Error: The specified directory for merging annotations (--distill_custom) does not contain any .tsv files: .`


There was a big bug in combine annotations - was missing some annotations but formatting was still okay.
I will need to re-run the test data for Maddie.... eventually




## Monday March 25
### Recap
- in a good place
- want to run some more "big" tests to make sure the current version is stable.
- Then, I want to merge the dev into main, to keep it, and then continue working on dev with more changes
  


### Testing conda profile:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-profile-testing-03252024 --threads 5 --call --rename --annotate --use_camper --use_dbcan --use_kegg --distill_topic carbon --distill_ecosystem 'eng_sys ag' -profile conda_slurm --slurm_node tardigrade
```
- Works!
	- added separate config files for each profile: singularity w/w/o slurm, conda w/wo slurm
		- set memory as a variable
		- for conda, all processes except quast(environment-quasy.yml) use the same environment.yml
	- memory set in nextflow.config
		- could add in further tuning for user
			- this could be the DRAM-lite mode..
	- running `-profile singularity_slurm` - now, server is busy though



### Tests to run over night:

76 HMP Set:
- need to edit command
```groovy
nextflow run -bg DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-13-dbs-CAD-03182024 --threads 8 --call --rename --annotate --use_uniref --use_kegg --use_merops --use_viral --use_pfam --use_camper --use_kofam --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -with-report -with-trace -with-timeline
```

162 Gary+Jared Set:
- need to edit command
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/gary_jared_test/ --outdir DRAM2-Gary-Jared-162-CAD-03152024 --threads 20 --call --rename --annotate --use_kegg --use_uniref --use_kofam --use_vog --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --slurm_node tardigrade -resume --taxa ../test_data/gary_jared_test/gary_jared_taxa_binqual_162.tsv --bin_quality ../test_data/gary_jared_test/gary_jared_taxa_binqual_162.tsv
```  



Creating logo:
```
- I need design logo for software. Software used for the distillation of microbial gene annotations. Theme design is like making alcohol - first getting the "raw" annotations and then we are "distilling" this information for the user. Thus, we want a microbe-themed and alcohol distilling. single beaker,distilling flask, laying in dirt. Microbial pathways, in rainbow liquid, inside flask dripping out the distill head liquid into separate color puddles.

    Designer](https://www.bing.com/images/create?FORM=GDPGLP "Designer")|1024 × 1024 jpg
 
    Content credentials
    
    Generated with AI ∙ March 25, 2024 at 7:18 PM
```



## Wednesday March 27
- Modified the product process a bit to accommodate Maddie's scripts
- Going to try running the conda version on Unity
	- Trying to debug why Unity is not working for both the [[Gene-Centric Pipeline]] and [[COMET]]


Unity testing:
```groovy
nextflow run DRAM2.nf --input_fasta test_data/2_subsampled/ --outdir DRAM2-HMP-2-profile-testing-03252024 --threads 5 --call --rename --annotate --use_camper --use_dbcan --use_kegg --distill_topic carbon --distill_ecosystem 'eng_sys ag' --taxa test_data/gtdbtk.bac120.summary.tsv --bin_quality test_data/checkm1-HMP-76-quality.tsv -profile conda_slurm --slurm_queue 'wrighton'
```
- Still cannot get Unity to work - Sent message to Unity staff


-----

## Monday April 1
- Made nextflow.config the only one (no nextflow-Server.config)
	- This was just generalizing the profile config files
	- Need to run a small test

Test command for update - conda_slurm:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-profile-testing-01012024 --threads 10 --call --rename --annotate --use_camper --use_dbcan --use_kegg --distill_topic carbon --distill_ecosystem 'eng_sys ag' -profile conda_slurm --slurm_node tardigrade
```
Results:
- bbmap=39.06 did not work and complained about samtools
	- removed version specification and trying again
- Finished! Worked!
	- 16 minutes
  
Test command for update - singularity_slurm:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/2_subsampled/ --outdir DRAM2-HMP-2-profile-testing-01012024 --threads 10 --call --rename --annotate --use_camper --use_dbcan --use_kegg --distill_topic carbon --distill_ecosystem 'eng_sys ag' -profile singularity_slurm --slurm_node tardigrade --slurm_queue 'wrighton-hi,wrighton-low'
```
Results:
- Worked! only 3 minutes!
	- Seems conda is quite slower...


Rebuilt the singularity image to have the dependencies for maddie:
`DRAM2-Nextflow-Main-Container-March262024-V5.sif`
- [ ] Need to test



## Thursday April 4
- Found out there is a queueSize parameter
	- Not sure yet exactly how it works but going to test on rename

Rename test command for queueSize:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-ranme-queueSize-test-04042024 --threads 1 --rename --slurm_node tardigrade --queue_size 1 --max_forks_single_cpu 20 -profile singularity_slurm
``` 
Results
- Who won, queueSize or maxForks?
	- queueSize won - only 1 rename occurred at a time
- How does this work when annotating a lot of stuff?
  
```groovy
nextflow run -bg DRAM2.nf --input_fasta ../test_data/HMP_DRAM1_test_data/ --outdir DRAM2-HMP-76-ranme-queueSize-test-04042024 --threads 10 --call --rename --annotate --use_merops --use_viral --use_pfam --use_camper --use_dbcan --use_methyl --use_canthyd --use_vog --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 1 --max_forks_single_cpu 20 -profile singularity_slurm --slurm_queue 'wrighton-hi,wrighton-low' --slurm_node tardigrade
```


-----


## Friday April 12

Working on TREES.

I think this genome should have the KOs:
`/home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/test_data/15_soil_genomes_DRAM1/single_dechloromonas/Dechloromonas_aromatica_RCB.fasta`


Testing comamnd:
```groovy
nextflow run -bg DRAM2.nf --input_fasta ../test_data/15_soil_genomes_DRAM1/single_dechloromonas/ --outdir DRAM2-Trees-test-d_aromatica-040122024 --threads 10 --call --rename --annotate --use_kofam --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 3 -profile singularity_slurm --slurm_queue 'wrighton-hi,wrighton-low' --slurm_node tardigrade --threads 10
```

- almost working!
  
  
  
Switching to a different test data set:
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/nar_nxr_dataset/ --outdir DRAM2-X-nar_nxr-trees-test-04122024 --threads 10 --call --rename --annotate --use_kegg --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 3 -profile singularity_slurm --slurm_queue 'wrighton-hi,wrighton-low' --slurm_node tardigrade -resume --trees nar_nxr
```
Results
- None - server was packed - moved to Riviera:

Switching Riviera:
- switching to kofam as well - copying kegg over now..
```groovy
nextflow run DRAM2.nf --input_fasta ../test_data/nar_nxr_dataset/ --outdir DRAM2-7-nar_nxr-trees-test-04122024 --threads 30 --call --rename --annotate --use_kofam --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 0 -profile conda_slurm -resume --trees nar_nxr --max_forks_single_cpu 10 --max_forks_hmm 10 --max_forks_user_cpu 10 --time '1h'
```



Gary Jared test dataset:
```groovy
nextflow run -bg DRAM2.nf --input_fasta ../test_data/gary_jared_dataset/ --outdir DRAM2-Gary-Jared-162-CAD-04122024 --threads 30 --call --rename --annotate --use_kegg --use_vog --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 0 -profile conda_slurm -resume --max_forks_single_cpu 10 --max_forks_hmm 10 --max_forks_user_cpu 10 --time '1h' --taxa ../test_data/gary_jared_dataset/gary_jared_taxa_162.tsv --bin_quality ../test_data/gary_jared_dataset/gary_jared_binqual_162.tsv -with-trace -with-timeline -with-report
```  


## Monday April 22
- Trees integration is going OK
	- Need actual nar nxr tree
		- Need all other trees in pplacer input format
	- Need to modify trees logic - when to run and what to run



15 Soil metagenomes for Kayla:
```groovy
nextflow run -bg DRAM2.nf --input_fasta ../test_data/15_soil_genomes_DRAM1/ --outdir DRAM2-15-Soil-MetaG-CAD-04222024 --threads 30 --call --rename --annotate --use_kegg --use_vog --use_pfam --use_merops --use_viral --use_camper --use_dbcan --use_methyl --use_canthyd --use_fegenie --use_sulfur --distill_topic default --distill_ecosystem 'eng_sys ag' --queue_size 0 -profile conda_slurm -resume --max_forks_single_cpu 10 --max_forks_hmm 10 --max_forks_user_cpu 10 --time '3h' -with-trace -with-timeline -with-report -profile conda_slurm
```  

--------

## Wednesday May 8
-  Need to work on trees despite having no actual correct tree data. I am hoping in 2 months I finally have a real tree.. we shall see. 
- I will use the tree Kayla said to use but it does not seem right at all.. 
	- This is from the paper: [Extraordinary phylogenetic diversity and metabolic versatility in aquifer sediment](https://www.nature.com/articles/ncomms3120#Sec17) supplementary dataset 5.

--------

## Thursday May 16

Generated fake input data for Maddie for Product/Viz:

Datasets:
`/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/Testing-Results/DRAM2-HMP-76-13-dbs-CAD-03122024`
`/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/Testing-Results/DRAM2-Gary-Jared-162-CAD-03182024`



Next, going to think about Strainer:
Purpose of Strainer:
- Pull out gene sequences based on user desires
- User can provide
	- single gene ids
	- EC numbers 
	- taxonomy
	- 

Just going to put this on pause as I do not have enough information to implement well.


--------


## Friday May 17

### Meeting with Kayla
- She is going to start writing up a mock manuscript

--------


## Monday May 20
Working on prepraing the databases.
I have `/internal/prepare_databases.py`
- Thus far I have camper and vogdb working 

------


## Tuesday May 21
- Got pfam integrated into the `prepare_databases.py` however, it is ~16GB so i am not going to run it now.. will pop on server later.
	- also, i think pfam has some special formatting?
- Merops working
- working on fegenie but the individual hmm files are giving me trouble.
  
------

## Wednesday May 22
- Working on cleaning the fegenie concatenated hmm file

Semi-working command for fegenie:
`python assets/internal/prepare_databases.py --output_dir ../DRAM2-preparedatabases-test-output --databases fegenie --threads 5 --verbose --download_date 05202024`

This seems to be correct but I want to run it again and double check.

This makes is 6/13 databases able to download and index. (7th would be Uniref but it is too big to test ATM)

------


## Wednesday May 29
 - Going to work on cant hyd and dbcan - DONE
 - Strange, I am not sure where the cant_hyd faa file came from -  they only have the hmms on the GitHub - added to [[DRAM2-NF Big Questions]]
 - Pfam testing at some point but it is V large.

------


## Thursday May 30
- Working on kofam
	- have separate hmms. will concatenate and then index

------


## Friday May 31
- kofam working

Ran test:
```bash
python assets/internal/prepare_databases.py --output_dir ../DRAM2-preparedatabases-test-output --databases cant_hyd dbcan fegenie merops methyl sulfur vogdb kofam --threads 5 --verbose --download_date 05202024
```
Results:
Worked great!!


Databases left:
- pfam - mmsqes (not sure where to get the mmseqs b/c a hmm download link is in the DRAM1 and DRAM2 code) some code written
- uniref90 - mmseqs (not even sure if we are doing this one yet)
- kegg - mmseqs (need to review formatting code and find Rich's downloaded latest kegg)
- viral - mmseqs



------


## Tuesday June 25
- Set up test runs to show kayla how DRAM2 works

Profiles:
```groovy
-profile singularity_slurm 
-profile conda_slurm 
-profile singularity
-profile conda
```

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




