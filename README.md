# Resilience-versus-Environmental-Stochasticity

This repository contains R code associated with the analyses outlined in **Cant et al. (2023). Recent exposure to environmental stochasticity does not determine the demographic resilience of natural populations**, published in Ecology Letters.

If you have any questions, don't hestitate to get in touch: *jic2@st-andrews.ac.uk*

## RScript details

***Demographic Data Extraction***
This script outlines the initial process of extracting the transient characterisitics and relevant metadata of matrix population models (MPMs) sourced from the COMPADRE and COMADRE demographic databases (https://compadre-db.org/). This script requires the contents of these databases to be provided as .RData files. This script also relies on functions contained within the script ***Transient functions***.

***Environmental data extraction***
Following demographic data extraction this script allows users to locate and extract records of the local precipitation and temperature regimes for all selected MPMs from the CHELSA climate database (https://chelsa-climate.org/). These abiotic time-series can then used to quantify measures of environmental stochasticity using the ***Quantifying Environmental Variability*** script.

***Phylogeneitic extraction*** & ***Tree tricker function***
These scripts enable the formulation of a population level phylogenetic tree comprising the taxonomic details of all selected MPMs.

***Phylogenetic Partial Least squares Analysis***
This script forms the main body of the analyses outlined in *Cant et al. (2023) Ecology Letters* and relies on the ***Phylogenetic Partial Least squares function*** script, to simultaneously assess the relationship between exposure to environmental stochasticity and demographic resilience, and the vital rate sensitivities of resistance, recovery and compensation.

The final two scripts ***PLS function for short-lived_long-lived analysis*** & ***PPLS analysis short vs Long-lived species*** are largely repititions of the two scripts outlined immediately above, but they explore how any findings are influenced by the relative longevity of the selected populations.
