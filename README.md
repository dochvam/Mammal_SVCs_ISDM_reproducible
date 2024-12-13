# README

This repository contains code for reproducing the analyses presented in 
"Mammal niches are not conserved over continental scales" by Goldstein et al. 
This README document provides an overview for how to run this code and also
for how to access and explore the published input data and results associated
with the manuscript. For full methodological details, see the main manuscript.

## Organization of code

This repository contains .R code for reproducing the model estimation pipeline
presented in the manuscript.

There are three core files in the main directory for reproducing modeling:
`run_models_lineage.R`, which executes the full model estimation pipeline for
species with phylogeographic information; `run_models_nolineage.R`, which
executes the full model estimation pipeline for all species, discounting
phylogeographic information; and `run_simulation.R`, which executes the
simulation-based validation exercise presented in the manuscript.
Two files are provided for summarizing model output: `summarize.R` executes
post-model-fitting summary analyses including evaluating support for the core
hypotheses and producing predictive maps of species density; and `paper_figs.R`
produces figures.

The file `necessary_packages.R` loads in all the necessary packages required for
all scripts. It is not necessary to execute this file but it is provided so you
can check that all necessary packages are loaded.

Be aware that analyses presented are very slow and took us weeks to execute in
full. In the model execution files, analyses are coded to execute in parallel
(across species or simulation scenario models) and the number of cores is
specified, though this number may not be appropriate for all systems.


Helper code is organized into the following directories:

- `setup_code/` contains general helper functions 
- `main_code_lineage/` contains code for executing the full model estimation workflow for 13 species associated with ancient genetic lineages
- `main_code_nolineage/` contains code for executing the full model estimation workflow for all 35 species studied, not including information about ancient genetic lineages
- `validation_sim/` contains helper code for the simulation exercise


## Retrieving auxiliary data

Files containing iNaturalist counts, camera trap detection histories,
and spatial covariate data are hosted on Zenodo as they exceed the 
allowable file size for Github. To retrieve these, go to the following link and download the files from Zenodo:

[https://doi.org/10.5281/zenodo.14219130](https://doi.org/10.5281/zenodo.14219130)

Place the unzipped files from the `model_inputs/` Zenodo subdirectory into the 
directory `intermediate/`, and the repository should then function as intended.

The auxiliary data also include a README.txt file that describe the contents of each file.
