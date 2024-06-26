# DRAM2 Trees development notes and background



## DRAM2 Trees background and motivation

The purpose of DRAM2 Trees is to reslove differences in genes which have identical annotations. This is not to say these annotations are incorrect, it is to say that the annotation databases are not able to resolve the annotations to the level a user may need to state certain things about the genes they have annotated. 

Example: KEGG has a single KO for the genes: narG, narZ, nxrA. This results in the user not knowing which of these genes is actually present within their samples. However, researchers have been able to build phylogenetic trees with sequences of known subunits related to nitrogen reductase. Thus, by placing these sequences, which have identical KEGG KOs, on a tree with know sequences, one can decipher which sequences is which which subunit. 

DRAM2 Trees is being developed to do exactly this: identify pre-determined gene identifiers, obtained through annotation via the various databases, and place the corresponding sequences of these called genes on a pre-build phylogenetic tree. Then, based on either a closest neighbor metric or a pre-determined distance metric, assign a correct annotation to this called gene and update the annotations TSV file.

DRAM2 Trees will have pre-built trees for the following types of genes: nar/nxr, dsr, mcra, pmoa/amoa, and hydrogenases.



## DRAM2 Trees development

### Current code status and how Trees is designed to be run

DRAM2 Trees should run if traits/adjectives is going to be run. This is due to the fact that the rules which assign traits/adjectives to input samples depends on the corrected annotations for the specified genes (nar/nxr, dsr, mcra, pmoa/amoa, and hydrogenases). However, depending on the computational demand, and time demand, to run these various trees, it may be beneficial to NOT run DRAM2 Trees in the event that user does not want to run traits/adjectives. Either way, the GitHub documentation should clearly state if the annotations the user gets have been corrected using phylogenetic trees.

The current state of DRAM2 Trees is only for the dmso tree and relies only on a nearest neighbor decision for calling the gene a given annotation. Then, DRAM2 parses the annotations TSV file and generates a new column named "tree_verified" 

### How to prepare the REFPKG

1) Obtain fasta sequences

2) use muscle to align the sequences:

`grep -c ">" dmso_fortree.fasta  ## there are 3830 seqs`

`muscle -in dmso_fortree.fasta -out dmso_fortree_aligned.fasta`

3) protpipeliner.py to run raxml and find optimal sub model

`protpipeliner.py -i dmso_fortree_aligned.fasta -t 10 -b 100 -m low -a T`

4) Determine the number of positions using seqkit:

`seqkit stats dmso_refs.fasta.al`

```
#Output:
#file                format  type     num_seqs  sum_len  min_len  avg_len  max_len
#dmso_refs.fasta.al  FASTA   Protein        86  160,046    1,861    1,861    1,861
```

5) Create the `.refpkg` using taxit:

`taxit create -P dmso_package -l 1861 -p dmso_refs.fasta.model -f RAxML_info.dmso_refs.fasta --stats-type RAxML`

### Future development notes