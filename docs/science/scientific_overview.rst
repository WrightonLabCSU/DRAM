DRAM2 Impetus
====================

Large microbiome data demands robust methods to determine the functional potential of the staggering diversity of bacteria, archaea, and viruses. For this purpose, we created `DRAM1 <https://github.com/WrightonLabCSU/DRAM>`_ (Distilled and Refined Annotations of Metabolism), a tool to analyze and curate genome annotations from multiple databases, summarizing the functions encoded by each genome. 

DRAM2 is a significant upgrade to DRAM1 that adapts to the larger scale of datasets and evolving complexity of research in the microbiome. DRAM2 provides more customization of the work-flow, including database selection and returned information in metabolism summaries. We addressed challenges in homology-based annotation by making use of phylogenetic trees and gene synteny, adding overall improvements to functional annotations. With these new features and a set of expertly-curated rules, DRAM2 can provide users with a per-genome list of metabolic traits or adjectives (e.g. Methanogen, Fermenter, Sulfate reducer) to allow for rapid genome profiling. 

Further, we have enabled SnakeMake pipeline integrations, allowing DRAM2 to run on a broader range of architectures with high performance. Like DRAM version 1, DRAM2 *will* be freely accessible in KBase, a point-and-click interface, alleviating the need for computational resources. This link to Kbase will also enable integration of DRAM2 annotations into other tools, such as genome scale metabolic models. We show that DRAM2 can surpass its predecessor by comparing both approaches on a large Wetland Microbiome database. This demonstrates the immediate value of DRAM2’s innovations and how they allow for a new understanding of microbiome function. 



