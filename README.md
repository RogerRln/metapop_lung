# Code for: Metapopulation model of phage therapy of an acute *Pseudomonas aeruginosa* lung infection

[![MATLAB Build](https://github.com/RaunakDey/metapop_lung/actions/workflows/blank.yml/badge.svg?event=push)](https://github.com/RaunakDey/metapop_lung/actions/workflows/blank.yml)

**Abstract**:
Infections caused by multi-drug resistant (MDR) pathogenic bacteria are a global health threat. Phage therapy, which uses phage to kill bacterial pathogens, is increasingly used to treat patients infected by MDR bacteria. However, the therapeutic outcome of phage therapy may be limited by the emergence of phage resistance during treatment and/or by physical constraints that impede phage-bacteria interactions in vivo. In this work, we evaluate the role of lung spatial structure on the efficacy of phage therapy for *Pseudomonas aeruginosa* infection. To do so, we developed a spatially structured metapopulation network model based on the geometry of the bronchial tree, and included the emergence of phage-resistant bacterial mutants and host innate immune responses. We model the ecological interactions between bacteria, phage, and the host innate immune system at the airway (node) level. The model predicts the synergistic elimination of a *P. aeruginosa* infection due to the combined effects of phage and neutrophils given sufficiently active immune states and suitable phage life history traits. Moreover, the metapopulation model simulations predict that local MDR pathogens are cleared faster at distal nodes of the bronchial tree. Notably, image analysis of lung tissue time series from wild-type and lymphocyte-depleted mice (n=13) revealed a concordant, statistically significant pattern: infection intensity cleared in the bottom before the top of the lungs. Overall, the combined use of simulations and image analysis of in vivo experiments further supports the use of phage therapy for treating acute lung infections caused by *P. aeruginosa* while highlighting potential limits to therapy given a spatially structured environment, such as impaired innate immune responses and low phage efficacy.

## Code usage

All the scripts used to generate main and supplementary figures are written in MATLAB (R2020b).

1. Run the script 'Generate_all_figs.m' inside the metapop_lung main folder to generate main and supplementary figures.
   1. In case you require more Java Heap Memory to run the 'Generate_all_figs.m' script, go to MATLAB -> Preferences -> General -> Java Heap memory, and increase the Java heap size, then click on apply and okay.
 
- [**metapop_code**](./metapop_code): directory contains all the scripts (.m) and functions necessary to generate metapopulation model simulations
- [**image_analysis**](./image_analysis): directory contains all the scripts (.m) and functions necessary to carry out imaging analysis of *P.a.* infected mice
- [**figures**](./figures): directory with saved figures after running [Generate_all_figs.m](./Generate_all_figs.m)
- [**data**](./data): directory with necessary data to run the metapopulation model scripts

**Note:**
I have pre-saved the data obtained after performing the robustness analysis of the metapopulation model inside the **data** folder. Hence, the scripts that generate the figures from the robustness analysis used the pre-saved data. To generate the data from scratch (it takes 13 hr per code), uncomment lines 45 to 73 of the script ['Fig5_heatmap_adsorption_rate_nlung.m'](./metapop_code/Fig5_heatmap_adsorption_rate_nlung.m) and lines 40 to 69 of the script ['FigS4_heatmap_mucin_nlung.m'](./metapop_code/FigS4_heatmap_mucin_nlung.m) both scripts are inside the [**metapop_code**](./metapop_code) folder.

