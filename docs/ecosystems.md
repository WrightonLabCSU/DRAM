# DRAM2 Ecosystems

DRAM2 Ecosystems provide ecosystem-specific summaries of microbial metabolic
potential. Each ecosystem organizes DRAM2 annotations into biologically
relevant processes and traits, allowing users to interpret microbial functions
within the context of a particular environment.

Rather than examining individual gene annotations alone, DRAM2 Ecosystems group
genes and pathways into higher-level metabolic traits relevant to specific
microbial systems. These outputs are designed to help users identify organisms
with the potential to contribute to important ecosystem processes and compare
functional potential across genomes, communities, samples, or environmental
conditions.

DRAM2 currently includes ecosystem-specific outputs for:

- **Biogeochemical Cycling** — environmental microbial communities and elemental cycling
- **Agricultural Systems** — managed soils and agricultural environments
- **Engineered Systems** — industrial, treatment, energy, and other managed microbial systems
- **Gut** — human and animal gastrointestinal microbiomes

The ecosystem selected should reflect both the environment being studied and
the biological questions being asked. Some metabolic processes occur in
multiple ecosystems, but each ecosystem organizes and interprets these
functions in a different biological context.

---

## Biogeochemical Cycling

### Ecosystem Overview

The Biogeochemical Cycling ecosystem is designed to characterize microbial
metabolic potential involved in the transformation of carbon, nitrogen,
sulfur, methane, hydrogen, metals, and other environmentally relevant
compounds. These processes connect microbial metabolism to the movement and
transformation of elements through environmental systems.

Environmental microbial communities frequently experience strong gradients
in oxygen availability, redox potential, organic carbon composition,
electron donors, and electron acceptors. Microorganisms respond to these
conditions through diverse respiratory, fermentative, carbon-fixation, and
energy-conservation strategies. The distribution of these metabolic
capabilities can therefore provide insight into the potential roles of
individual organisms within biogeochemical cycles.

The Biogeochemical Cycling ecosystem organizes annotations into major
environmental metabolic processes, including methane metabolism, carbon
fixation, carbon degradation, fermentation, hydrocarbon metabolism, sulfur
cycling, iron metabolism, hydrogen metabolism, photosynthesis, carbon monoxide
oxidation, electron transport, and environmental adaptations.

Use the **Biogeochemical Cycling** ecosystem when interpreting microbial
genomes or environmental -omics datasets where the primary goal is to
understand microbial contributions to elemental cycling and environmental
metabolism rather than processes specific to a host, agricultural system, or
engineered environment.

### Databases Required

- KEGG
- dbCAN
- MEROPS
- CAMPER
- FeGenie
- Sulfur
- Methyl
- Metals
- CARD
- dram_db
- TCDB

### Major Metabolic Processes

| Major Process | Included Subtopics |
| --- | --- |
| **Methane Metabolism** | Methanogenesis, methanogenic substrates, anaerobic methane oxidation, aerobic methane oxidation |
| **Carbon Fixation** | Calvin cycle, Wood-Ljungdahl pathway, reductive TCA cycle, 3-hydroxypropionate pathways, dicarboxylate/4-hydroxybutyrate pathways, reductive glycine pathway |
| **Carbon Degradation** | Polymer degradation, aromatic compound degradation, monomeric sugar utilization |
| **Fermentation** | Short-chain fatty acid production, alcohol production |
| **Hydrocarbon Metabolism** | Aerobic and anaerobic hydrocarbon degradation |
| **Hydrogen Metabolism** | Hydrogenases and hydrogen-dependent energy metabolism |
| **Sulfur Cycling** | Sulfate and sulfur reduction, sulfur oxidation, other sulfur transformations |
| **Iron Metabolism** | Iron oxidation, iron reduction |
| **Photosynthesis** | Photosynthetic metabolic potential |
| **Carbon Monoxide Oxidation** | CO oxidation |
| **Electron Transport** | Components of the electron transport chain and energy conservation |
| **Environmental Adaptations** | Environmentally relevant physiological and stress-response traits |

For the exact genes, identifiers, and rules used to define each metabolic
trait, consult the Biogeochemical Cycling ecosystem ruleset.

### Interpretation Guidance

DRAM2 Biogeochemical Cycling outputs summarize genomic evidence for metabolic
processes involved in environmental element and energy cycling. Positive trait
calls indicate that the genomic evidence specified by the corresponding rule
was detected. Because some metabolic processes can be represented by marker
genes while others require multiple enzymes, complexes, or pathway
components, users should consult the underlying rule and gene annotations when
interpreting individual calls.

Pathway completeness provides information about the components of a defined
metabolic pathway that were recovered. More complete pathways generally
provide stronger evidence for metabolic potential, while partial pathways
should be evaluated in the context of genome completeness and the biological
requirements of the pathway.

Trait predictions represent **metabolic potential rather than activity**.
For example, detection of genes associated with methanogenesis, sulfur
oxidation, carbon fixation, or hydrocarbon degradation indicates that an
organism may have the genetic capacity for that process. It does not establish
that the process is occurring under the sampled environmental conditions.

Biogeochemical interpretations should therefore consider environmental
context, particularly redox conditions, oxygen availability, substrate
availability, and the presence of appropriate electron donors and acceptors.
Multiple organisms may also contribute different steps of the same
biogeochemical transformation, so community-level interpretation should not
assume that every pathway must occur within a single genome.

---

## Agricultural Systems

### Ecosystem Overview

Agricultural ecosystems are managed environments in which plant productivity,
livestock, soil properties, climate, and management practices interact to
shape microbial communities and their metabolic potential. Across cropping
and grazing systems, agricultural soils are highly heterogeneous and
experience repeated changes in resource availability, redox conditions,
moisture, disturbance, and plant inputs.

Microbial metabolism regulates many processes that determine agricultural
productivity and environmental outcomes. Microorganisms transform plant
residues and soil organic matter, cycle nitrogen, phosphorus, and sulfur,
mediate greenhouse gas production and consumption, and contribute to nutrient
acquisition and stress responses. These functions influence soil carbon
persistence, nutrient availability, fertilizer-use efficiency, plant
productivity, and the potential for nutrient loss.

The Agricultural Systems ecosystem organizes annotations into processes
relevant to managed soils, including carbon processing, inorganic and organic
nitrogen transformations, nutrient uptake and transport, phosphorus
mobilization, sulfur metabolism, phytohormone-related functions, stress
tolerance, and potential pesticide degradation. It also connects genomic
annotations with commonly measured soil enzyme activities.

Use the **Agricultural Systems** ecosystem when interpreting metagenomic
assemblies or genomes from managed soils and other agricultural environments,
particularly when questions related to nutrient cycling, carbon
transformation, plant-associated functions, or microbial responses to
management are central to the study.

### Databases Required

- KEGG
- dbCAN
- MEROPS
- CAMPER
- FeGenie
- Sulfur
- Methyl
- Metals
- CARD
- dram_db
- TCDB

### Major Metabolic Processes

| Major Process | Included Subtopics |
| --- | --- |
| **Carbon Processing** | Cellulose degradation, hemicellulose degradation, lignin and aromatic compound oxidation, starch degradation, pectin degradation, sugar utilization, short-chain fatty acid metabolism, alcohol metabolism, carbon dioxide fixation, methane oxidation, multicopper oxidases, laccases, lignin peroxidases |
| **Inorganic Nitrogen Transformations** | Urea degradation, nitrogen fixation, nitrification, aerobic ammonia or ammonium oxidation, nitrite oxidation, anammox, nitrate reduction, DNRA, denitrification |
| **Organic Nitrogen Transformations** | Extracellular peptidases, aspartic peptidases, cysteine peptidases, metallopeptidases, asparagine peptidases, glutamic peptidases, serine peptidases, threonine peptidases, chitin degradation, amidases, formamidases, ureases |
| **Nitrogen Uptake and Transport** | Nitrate and nitrite transport, ammonium transport, urea transport |
| **Phosphorus Transformations** | Inorganic phosphorus solubilization, phosphatase-mediated organic phosphorus mineralization, phytase-mediated phosphorus mineralization, phosphonate mineralization |
| **Phosphorus Uptake and Transport** | Organic phosphorus transport, inorganic phosphate transport |
| **Sulfur Metabolism** | Inorganic and organic sulfur transformations, sulfur oxidation and reduction, sulfur assimilation, sulfur compound transport |
| **Enzyme Assay Analogs** | Beta-D-glucosidase, beta-D-xylosidase, arylsulfatase, leucine aminopeptidase, N-acetyl-beta-D-glucosaminidase, phosphatase, cellobiohydrolase, alpha-1,4-glucosidase, peroxidase |
| **Phytohormone-Related Functions** | Auxin/IAA biosynthesis and inactivation, cytokinin biosynthesis and degradation, ACC deaminase, salicylic acid biosynthesis and degradation |
| **Stress Tolerance and Persistence** | Compatible solutes, oxidative-stress responses, xenobiotic detoxification, heavy-metal transformation and resistance, EPS production, biofilm production |
| **Potential Pesticide Degradation** | Glyphosate degradation, pyrethroid degradation, organophosphate degradation, atrazine degradation |

For the exact genes, identifiers, and rules used to define each metabolic
trait, consult the Agricultural Systems ecosystem ruleset.

### Interpretation Guidance

DRAM2 Agricultural Systems outputs summarize genes and gene combinations
associated with agriculturally relevant microbial functions. A positive result
indicates that the genomic evidence required by the corresponding rule was
detected. Some functions are represented by individual markers, whereas
others require combinations of enzymes or pathway components. Users can
examine the underlying annotations and rule structure to determine the
evidence supporting each call.

Pathway completeness describes the proportion of required pathway components
detected in an assembly or genome. More complete pathways generally provide
stronger evidence for the corresponding metabolic capacity. Partial pathways
identify which components were recovered but should be interpreted alongside
genome quality because assembly and binning gaps can result in apparent gene
absence.

Trait predictions describe **genetic potential rather than demonstrated
activity**. For example, a genome containing carbohydrate-active enzymes and
extracellular degradation systems may have the potential to participate in
plant-residue decomposition. Phosphate-solubilization or
organic-phosphorus-mineralization markers may indicate potential involvement
in phosphorus mobilization, while compatible-solute, oxidative-stress, or
metal-resistance traits may indicate adaptations to environmental stress.

Comparing these functional profiles across genomes, communities, management
systems, or environmental gradients can help connect microbial metabolic
potential to broader patterns in carbon cycling, nutrient availability, plant
interactions, and soil responses to agricultural management.

---

## Engineered Systems

### Ecosystem Overview

Engineered systems are human-designed or human-managed environments such as
energy extraction operations, including hydrocarbon and geothermal systems,
wastewater treatment plants, bioreactors and anaerobic digesters,
contaminated-site remediation, industrial pipelines, and other built
infrastructure. Unlike natural ecosystems, engineered systems are defined by
human intervention and often have controlled or manipulated physical and
chemical parameters such as redox conditions, temperature, and substrate
loading.

Microorganisms can contribute to both undesirable outcomes, such as corrosion
and biofilm formation, and desirable outcomes, such as contaminant removal,
substrate transformation, and product generation. Microbial communities in
engineered systems are therefore frequently shaped by intentional system
management in addition to environmental gradients and nutrient availability.

The major processes captured by the Engineered Systems ecosystem include
sulfur cycling, methane production and consumption, hydrogen metabolism,
nitrogen cycling, metal cycling and resistance, fermentation, carbon and
hydrocarbon degradation, biofilm formation, and environmental stress
responses.

Microbial metabolism is central to engineered systems because the collective
biochemistry of the microbial community can directly affect system
performance. Microorganisms can degrade organic or hydrocarbon substrates,
cycle nutrients, produce or consume methane and hydrogen, transform metals,
and form or resist biofilms. Changes in these processes can directly
influence system efficiency, stability, product formation, contaminant
removal, or infrastructure integrity.

Use the **Engineered Systems** ecosystem when analyzing -omics data from
managed or industrial microbial systems where the goal is to link microbial
metabolic potential to system performance, troubleshooting, or optimization.

### Databases Required

- KEGG
- Cant-Hyd
- dram_db
- FeGenie
- Sulfur
- dbCAN
- TCDB

### Major Metabolic Processes

| Major Process | Included Subtopics |
| --- | --- |
| **Sulfur Cycling** | Sulfide production, sulfur oxidation, sulfur transformations |
| **Methane Cycling** | Methane production, methane consumption |
| **Hydrogen Use** | FeFe hydrogenases, NiFe hydrogenases |
| **Nitrogen Cycling** | Nitrification, anammox, comammox, denitrification, DNRA |
| **Metals** | Metal cycling, iron acquisition, metal detoxification |
| **Fermentation End Products** | Organic acid production including acetate, propionate, butyrate, lactate, formate, and succinate; alcohol production including ethanol and butanol |
| **Carbon Degradation** | Aerobic hydrocarbons, anaerobic hydrocarbons, polymers, simple sugars, glycerol, ethylene glycol, methanol |
| **Biofilms** | EPS production, quorum sensing, motility |
| **Environmental Stress Response** | Osmoregulation, heat and cold stress, antimicrobial resistance, biocide resistance |

For the exact genes, identifiers, and rules used to define each metabolic
trait, consult the Engineered Systems ecosystem ruleset.

### Interpretation Guidance

DRAM2 Engineered Systems outputs summarize genomic evidence for metabolic
functions relevant to managed and industrial microbial systems. In the
heatmap, presence indicates that a genome or MAG contains the genes or
subunits required by the rule defining that metabolic trait. Rules differ
among traits, so users should consult the Engineered Systems ruleset to
determine the evidence required for an individual call. The complete
Engineered Systems summary file can be used to examine the individual genes
and subunits underlying these trait-level summaries.

Pathway completeness describes how much of the defined gene set for a pathway
was recovered. More complete pathways generally provide stronger genomic
support for a metabolic capability, while partial pathways should be
interpreted alongside genome or MAG completeness because assembly and binning
gaps can result in missing genes.

Trait predictions represent **metabolic potential rather than activity**.
For example, several genomes may encode the potential for acetate production
without acetate accumulating in the system. Acetate may be consumed by other
microorganisms, influenced by abiotic processes, or generated through pathways
that also participate in central metabolism.

Interpretation should therefore consider system conditions such as redox
state, substrate availability, temperature, pH, and biotic and abiotic
interactions. Where available, genomic predictions should be integrated with
metatranscriptomic, metabolomic, chemical, and process measurements to connect
microbial metabolic potential with engineered-system performance.

---

## Gut

### Ecosystem Overview

Gut ecosystems are internal, host-regulated environments in which host diet,
physiology, immune activity, and microbial community interactions shape
microbial metabolic potential. Gut environments are compartmentalized and
range from anaerobic to microaerophilic conditions. Microorganisms encounter
dietary substrates, host-derived compounds, host tissues, immune defenses,
and metabolites produced by other members of the microbial community.

The major processes captured by the Gut ecosystem include dietary carbon
processing, fermentation, redox and alternative electron acceptor use,
colonization and persistence, host interactions, and pathogenesis and
virulence.

Microbial metabolism in the gut can directly influence host physiology.
Fermentation of dietary and host-derived carbon generates short-chain fatty
acids and other metabolites that affect host energy harvest and gut barrier
integrity. Bile acid transformation, vitamin biosynthesis, and other
host-interaction pathways can influence host physiology and immune signaling.
Colonization traits influence whether microorganisms persist against host
clearance and competing taxa, while pathogenesis and virulence potential
provide information about possible adverse host-microbe interactions.

Use the **Gut** ecosystem when interpreting metagenomic or
metatranscriptomic data from human or animal gastrointestinal systems. The
human gastrointestinal tract is the primary point of reference, but this
ecosystem can also be used for rumen and other animal systems. It is most
useful for comparing functional capabilities across genomes, identifying
host-relevant traits, and summarizing host-microbe and microbe-microbe
interactions.

### Databases Required

- KEGG
- dbCAN
- antiSMASH
- FeGenie
- Metals
- Sulfur
- MEROPS
- dram_db

### Major Metabolic Processes

| Major Process | Included Subtopics |
| --- | --- |
| **Dietary Carbon Processing** | Resistant starch degradation, fiber degradation, glycan and mucin degradation, sugar utilization |
| **Fermentation** | SCFA production, lactate, succinate, ethanol, hydrogen metabolism, methanogenesis, acetogenesis |
| **Redox** | Oxygen utilization, nitrate reduction, sulfate reduction, TMAO reduction, fumarate respiration |
| **Colonization and Persistence** | Biofilms, EPS production, motility, quorum sensing, stress tolerance |
| **Host Interaction** | Bile acid metabolism, vitamin biosynthesis, polyphenol metabolism, drug metabolism, TMA production, TMA demethylation |
| **Pathogenesis and Virulence** | Toxins, adhesion and invasion, secretion systems, iron acquisition, hemolysins, immune evasion, oxidative-stress defense, biofilm-associated virulence |

For the exact genes, identifiers, and rules used to define each metabolic
trait, consult the Gut ecosystem ruleset.

### Interpretation Guidance

DRAM2 Gut outputs summarize potential metabolic functions based on genes and
gene combinations detected within a genome or MAG. Some rules use individual
marker genes, whereas others represent multi-step pathways. Users should
examine the rule structure and supporting annotations when evaluating the
strength of an individual metabolic call.

Pathway completeness reflects the number of components required by a rule
that are detected in the genome. It does not necessarily indicate that a
pathway is biologically complete. Incomplete calls should also be interpreted
alongside genome or MAG quality because assembly and binning gaps may result
in missing genes.

Trait predictions describe **genetic potential and do not demonstrate that a
function is expressed under a given condition**. For example, a genome
containing fiber- or mucin-degradation genes together with SCFA-production
markers may have the potential to break down complex carbon and produce
fermentation products. Bile acid or TMA-related markers indicate potential
involvement in transformation of host-associated metabolites. Motility,
biofilm, and stress-tolerance traits can provide information about the
potential for colonization and persistence.

Virulence-related traits should also be interpreted as genomic potential
rather than evidence that an organism is pathogenic. Individual virulence
markers should be evaluated in the context of the complete genome, taxonomy,
other host-interaction traits, and available experimental or clinical
information.

Comparing functional profiles across genomes or communities can help connect
microbial metabolic potential to broader patterns in dietary substrate use,
microbial interactions, host-associated metabolism, colonization, and host
physiology.
