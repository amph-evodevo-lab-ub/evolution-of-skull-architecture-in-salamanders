Title of the study: Evolution of cranial architecture and inter-element articulations in salamanders of the family Salamandridae

Summary
Mobility of skull elements relative to the braincase is a widespread vertebrate trait with major consequences for feeding performance, though its macroevolutionary dynamics remain poorly understood in amphibians. Using high-resolution micro-CT data from 278 specimens representing 54 species across 19 salamandrid genera (family Salamandridae), we quantified variation in ten cranial articulations over four skull regions and constructed a composite cranial articulation index. The articulation index varied substantially across species, ranging from largely separated (open) to fully consolidated (fused) articulations. Comparative analyses revealed that articulations of the snout and the skull roof evolved as an integrated module, whereas palatal and suspensorial articulations followed independent evolutionary trajectories, indicating a mosaic pattern of evolution in cranial articulation. Most articulations and the composite articulation index were best explained by Ornstein–Uhlenbeck models, consistent with evolution toward constrained or lineage-specific optima. We tested relationships between cranial articulation and skull size, skull robustness (bone thickness), and ecological preferences (length of the aquatic period). Contrary to expectations, our results indicate that cranial articulation in Salamandridae evolved as a modular system independent of these traits. These findings suggest that salamandrid cranial architecture is maintained near a structural optimum within a wider spectrum presumably governed by the complex synergy of developmental, phylogenetic and ecological constraints.

This README file include explanation of data and R code (R v. 4.5.2) used for analyses of evolution of cranial architecture and inter-element articulations in family Salamandridae using high-resolution micro-CT data from 278 specimens representing 54 species across 19 salamandrid genera.

Input Files
* R_data_salamandridae.csv   Scores of inter-bone articulations, cranial articulations index (CAI), skull robustness (corrected_thickness), skull size (CS) and ecological preferences (aquatic_period) averaged per species  which were used for further analyses
* Tree_salamandridae_cranial_articulations.nex   Time-calibrated molecular phylogeny from Stewart and Wiens (2025) including a resolved branching order for the Eurasian crested newts 

Dependencies
Install the following R libraries before running analyses:
library (ape) (version used for analyses v. 5.0)
library (phytools) (v. 2.0)
library(geiger) (v. 2.0)
library(dplyr) (v. 1.2.1)
library(igraph) (v. 2.2.1)
library(scales) (v. 1.4.0)Workflow Overview

Data Import & Preprocessing
* Read phylogenetic tree
* Read character states and mean values for each species
* Prune phylogenetic tree only to contain species with character states data

Map character states onto phylogeny with ancestral states
* Prepare results table
* Make a loop to analyze all traits: calculate ancestral states and map each character separately
Output:
* Table of ancestral state values for each character under Brownian motion (BM) and Ornstein Uhlenbeck (OU)
* Graph for each character with states mapped onto phylogeny along with OU (best supported) ancestral state  (presented in manuscript as Figure 3, supplementary figures S1-S13)

Calculate phylogenetic signal for each character
* Prepare results table
* Make a loop to analyze all traits: calculate OU ancestral states foe each character separately 
* False discovery rate (FDR) correction across all p-values
Output:
* Table of phylogenetic signal values for each trait quantified with Pagel's ? and Blomberg's K, p-values and FDR corrected p-values (presented in manuscript as Table 1)

Evolutionary model check
* Prepare results table
* Make a loop to analyze all traits fitting evolutionary models (BM, OU and early-burst EB) to phylogenetic trees for each character separately and evaluate model support using Akaike Information Criterion (AICc) and Akaike weights (w)
Output:
* Table of ancestral (root) states and support for all three models as AICc and Akaike weight for each character (presented in manuscript as Table 2)

Phylogenetically independent contrasts (PICs)
* Prepare results table
* Make a loop to analyze PICs for all trait pairs with statistical significance corrected for multiple comparisons using Benjamini-Hochberg procedure (BH) 
* Visualisation for cranial inter-bone articulations  
Output:
* Table of PIC regressions with the observed strength of association (R2), p-values and BH-corrected p-values among individual inter-bone articulations, cranial articulations index (CAI), skull size (CS) and skull robustness (SR) (presented in manuscript as Supplementary Table 1)
* Figure of network of statistically significant correlations among examined cranial articulation traits based on PICs where connecting edges represent R2 for statistically significant associations 

Citation
If using this workflow in your publication, please cite:
[suppressed for anonymity will be added upon acceptance]

Contact
For issues within script:
[suppressed for anonymity will be added upon acceptance]
For data requests:
[suppressed for anonymity will be added upon acceptance]

