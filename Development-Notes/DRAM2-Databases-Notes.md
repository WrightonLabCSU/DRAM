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