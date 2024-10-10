# DRAM2 Databases

## Database background

DRAM1 required the users to prepare their own databases. While this is a good idea, it requires the user to have >300GB memory to do so. This requirement immediately reduces the amount of potential DRAM users down significantly. Thus, to avoid this, DRAM2 aims to provide pre-built and formatted databases to users via GLOBUS. 

## DRAM2 Description database

DRAM1 and DRAM2 both rely on a pre-built descriptions SQL database which holds gene annotation metadata. This database is used, after annotation, to attach metadata to the annotation TSV for each called gene which was annotated. As with the DRAM2 databases below, this must be re-built or recycled from the DRAM1 code and must be rebuilt in the same way as DRAM2 currently functions with the DRAM1 build descriptions database. 

The descriptions database resides in:
`./databases/db_descriptions'

The descriptions database location is linked within the `nextflow.config` on line 195:

```
        // SQL annotation descriptions database
            sql_descriptions_db = "./databases/db_descriptions/description_db.sqlite"
```

## Current DRAM2 Database implementation

Unfortunately DRAM2 database implementation is minimal and thus will give the developer some options. Currently, within `./assets/internal/` there is a start to download and format the databases. While this may be a good path forward, one might also choose to re-use the DRAM2-archived code for `prepare_databases` or use the DRAM1 code. 

## DRAM2 database functionality

In the end, databases need to be formatted as they were in DRAM1. This requirement is due to the fact that, upon DRAM2 development, DRAM1 formatted databases, and the descriptions database, were relied upon. Thus, the format of these must be recapitulated in order for DRAM2 to function as it is currently functioning. 

## KEGG side-note

Rich was able to download the lasted version of KEGG before the Wrighton Lab's license expired (~April 2024). However, this does not include the pre-formatting which is needed to generate the DRAM2 databases. Thus, one must learn how DRAM1 formats the KEGG database into a usable format for input into DRAM1 `prepare_databases`. This is a large unknown in DRAM2 development. 


## Database desired structure and naming scheme
```bash
(base) zenith /home/projects-wrighton-2/Pipeline_Development/DRAM2-Nextflow/DRAM2-NF/databases: ls */*
dbcan/dbcan.fam-activities.tsv         kegg/kegg.mmsdb.idx           merops/merops.mmsdb_h.index     uniref/uniref.mmsdb             viral/viral.mmsdb.idx
dbcan/dbcan.fam.subfam.ec.tsv          kegg/kegg.mmsdb.idx.dbtype    merops/merops.mmsdb.idx         uniref/uniref.mmsdb.dbtype      viral/viral.mmsdb.idx.dbtype
dbcan/dbcan.h3f                        kegg/kegg.mmsdb.idx.index     merops/merops.mmsdb.idx.dbtype  uniref/uniref.mmsdb_h           viral/viral.mmsdb.idx.index
dbcan/dbcan.h3i                        kegg/kegg.mmsdb.index         merops/merops.mmsdb.idx.index   uniref/uniref.mmsdb_h.dbtype    viral/viral.mmsdb.index
dbcan/dbcan.h3m                        kegg/kegg.mmsdb.lookup        merops/merops.mmsdb.index       uniref/uniref.mmsdb_h.index     viral/viral.mmsdb.lookup
dbcan/dbcan.h3p                        kegg/kegg.mmsdb.source        merops/merops.mmsdb.lookup      uniref/uniref.mmsdb.idx         viral/viral.mmsdb.source
dbcan/dbcan.hmm                        kofam/kofam_ko_list.tsv       merops/merops.mmsdb.source      uniref/uniref.mmsdb.idx.dbtype  vogdb/vog_annotations_latest.tsv.gz
db_descriptions/description_db.sqlite  kofam/kofam_profiles.hmm      methyl/methyl.mmsdb             uniref/uniref.mmsdb.idx.index   vogdb/vog_latest_hmms.hmm
fegenie/fegenie.hmm                    kofam/kofam_profiles.hmm.h3f  methyl/methyl.mmsdb.dbtype      uniref/uniref.mmsdb.index       vogdb/vog_latest_hmms.hmm.h3f
fegenie/fegenie_iron_cut_offs.txt      kofam/kofam_profiles.hmm.h3i  methyl/methyl.mmsdb_h           uniref/uniref.mmsdb.lookup      vogdb/vog_latest_hmms.hmm.h3i
fegenie/old_entry_names.txt            kofam/kofam_profiles.hmm.h3m  methyl/methyl.mmsdb_h.dbtype    uniref/uniref.mmsdb.source      vogdb/vog_latest_hmms.hmm.h3m
kegg/kegg.mmsdb                        kofam/kofam_profiles.hmm.h3p  methyl/methyl.mmsdb_h.index     viral/viral.mmsdb               vogdb/vog_latest_hmms.hmm.h3p
kegg/kegg.mmsdb.dbtype                 merops/merops.mmsdb           methyl/methyl.mmsdb.index       viral/viral.mmsdb.dbtype
kegg/kegg.mmsdb_h                      merops/merops.mmsdb.dbtype    methyl/methyl.mmsdb.lookup      viral/viral.mmsdb_h
kegg/kegg.mmsdb_h.dbtype               merops/merops.mmsdb_h         methyl/methyl.mmsdb.source      viral/viral.mmsdb_h.dbtype
kegg/kegg.mmsdb_h.index                merops/merops.mmsdb_h.dbtype  sulfur/sulfur.hmm               viral/viral.mmsdb_h.index

camper/hmm:
camper.hmm  camper.hmm.h3f  camper.hmm.h3i  camper.hmm.h3m  camper.hmm.h3p  camper_hmm_scores.tsv

camper/mmseqs:
camper.mmsdb         camper.mmsdb_h         camper.mmsdb_h.index  camper.mmsdb.idx.dbtype  camper.mmsdb.index   camper.mmsdb.source
camper.mmsdb.dbtype  camper.mmsdb_h.dbtype  camper.mmsdb.idx      camper.mmsdb.idx.index   camper.mmsdb.lookup  camper_scores.tsv

canthyd/hmm:
cant_hyd.h3i  cant_hyd.hmm  cant_hyd.hmm.h3f  cant_hyd.hmm.h3i  cant_hyd.hmm.h3m  cant_hyd.hmm.h3p  cant_hyd_HMM_scores.tsv  CANT_HYD_HMM_scores.tsv

canthyd/mmseqs:
cant_hyd_BLAST_scores.tsv  cant_hyd.mmsdb.dbtype  cant_hyd.mmsdb_h.dbtype  cant_hyd.mmsdb.idx         cant_hyd.mmsdb.idx.index  cant_hyd.mmsdb.lookup
cant_hyd.mmsdb             cant_hyd.mmsdb_h       cant_hyd.mmsdb_h.index   cant_hyd.mmsdb.idx.dbtype  cant_hyd.mmsdb.index      cant_hyd.mmsdb.source

pfam/hmm:
Pfam-A.hmm.dat

pfam/mmseqs:
pfam.mmsmsa.dbtype  pfam.mmspro         pfam.mmspro_h         pfam.mmspro_h.index  pfam.mmspro.idx.dbtype  pfam.mmspro.index
pfam.mmsmsa.index   pfam.mmspro.dbtype  pfam.mmspro_h.dbtype  pfam.mmspro.idx      pfam.mmspro.idx.index

```