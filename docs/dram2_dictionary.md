# DRAM2 Dictionary

**ANNOTATE**: a DRAM2 function that pulls information from multiple databases because each database has different strengths, weaknesses, and purposes. Users are provided with results from all databases in the `raw-annotations.tsv` file, allowing them to use their knowledge of the system, organisms, and research questions to determine the most appropriate annotation for their genomes. This flexibility makes DRAM2 highly customizable to individual research needs. For a complete description of ANNOTATE see here: 

**SUMMARIZE**: a DRAM2 function that creates a summary of gene annotations by Topic and Ecosystem. DRAM2 organizes annotations into curated categories that highlight genes most relevant to specific metabolic functions, environments, and research applications. For a complete description of the SUMMARIZE function see here:

**QC**: Quality Control- QC follows MIMAG (*Minimum Information about a Metagenome-Assembled Genome*) standards, which define the minimum information required for reporting and publishing MAGs.


**Ecosystem** : describes a curated set of genes associated with pathways and processes relevant to a major area of microbial research. These include the following:
- **Agriculture (Ag)**: Genes associated with agriculturally relevant functions, including:
  - Nitrogen fixation and mineralization
  - Soil organic matter (SOM) degradation
  - Phytohormone production
  - Other functions that support microbial life in the rhizosphere
- **Engineered Systems**: Genes related to microbiological processes occurring in human-made environments, such as:
  - Waste treatment systems
  - Industrial bioprocesses
  - Fracking wells
- **Gut**: Genes associated with microbial functions in mammalian gut microbiomes.
- **Biogeochemical Cycling**: Genes involved in major elemental cycles, including:
  - Carbon cycling
  - Nitrogen cycling
  - Sulfur cycling
  - Phosphorus cycling
- **Water**: Genes associated with microbial processes in aquatic environments.


**Topic**: Topics are broad categories of metabolic activity used to organize information in the DRAM2 SUMMARIZE output. These include the following: 
- **Energy Acquisition / Bioenergetics** Describes how an organism obtains energy, including genes involved in:
  - Aerobic respiration
  - Anaerobic respiration
  - Other energy-generating pathways
- **Assimilation & Cofactor Metabolism** Describes which substrates an organism can acquire and utilize, including:
  - Nutrient uptake systems
  - Transporters
  - Substrate cleavage mechanisms
  - Cofactor biosynthesis and utilization
- **Cellular Machinery** Describes genes involved in core cellular processes, including:
  - Motility
  - Cell division
  - Housekeeping functions
  - Other fundamental cellular activities
- **Environmental Interaction & Adaptation** Describes how an organism interacts with and responds to its environment, including genes related to:
  - Stress tolerance
  - Quorum sensing
  - Biofilm formation
  - Environmental adaptation mechanisms


**Trait** : broad, genome-level metabolic function assigned to a MAG based on the presence of curated gene sets. Traits provide a high-level overview of an organism's potential ecological roles and metabolic capabilities, and include some of the following:
  - Aerobic: ≥50% of NADH dehydrogenase and at least one subunit of low affinity oxidase
  - Microaerophilic: ≥50% of NADH dehydrogenase and at least one subunit of high affinity oxidase
  - N fixer: nifHDK, anfG, or vnfG
     
For more information on DRAM2 **Traits**  see here:


**VIZUALIZE** : a DRAM2 function that generates an interactive heatmap for each the presence of genes associated with each **Trait**

*This is a living document, and is subject to modification by DRAM2 authors. Please reach out if there are any definitions you would like to see in this document!*
