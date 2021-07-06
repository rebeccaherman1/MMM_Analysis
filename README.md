# MMM_Analysis
Multi-Model Mean Analysis of Climate Simulations for Detection and Attribution

# Data

## Navigating Data Files
All observations and model output can be found in the [data](/data) folder, sorted by variable. Variables used for this analaysis include [pr](/data/pr) and [ts](/data/ts/MMM_data). All variables are averaged seasonally over JAS and area-averaged over the Sahel or some other ocean basin. 
Observation files are titled "observations.mat" while simulation files end with "\_all.mat". 
Observations are categorized by source, while simulations are categorized by institution, model, and run, and titled by the type of simulation (listed below).
* historical experiments (radiative forcings)
    - h: CMIP5 historical "ALL"  
    - a: CMIP5 historical Anthropogenic Aerosols "AA"  
    - n: CMIP5 historical Natural Forcings "NAT"  
    - g: CMIP5 historical Greenhouse Gases "GHG"  
    - cmip6_: same experiments, but in CMIP6.
* amip-piF: CMIP6 AMIP simulations forced with SST and preindustrial radiative forcings
* amip-hist: CMIP6 AMIP simulations forced with SST and historical ALL radiative forcings, which begin by 1850  

## Data Scripts
Below are the scripts used to generate these data files. 

* Observations and CMIP5 (IRIDL):
    - [Sahel_1_save_data_ingrid.m](/Sahel_1_save_data_ingrid.m): used to download observations to one file and CMIP5 simulations to another.
    - [institutions_cmip5.mat](/data/institutions_cmip5.mat): hand-made file used to identify CMIP5 institutions. Contains model name abbreviations (MODELS) and their corresponding "umbrella" research institution abbreviations (ABBREV) and full names (INSTITUTION) for all models used in Herman et al 2020 and followup work.
    - Model Lists: Files ending in "\_models.txt" contain lists of the files iterated over when downloading data using Sahel_1_save_data. Some of these files are currently missing from this repo.
* CMIP6 (cloud): 
    - [Download2NetCDF.ipynb](/Download2NetCDF.ipynb): used to download CMIP6 simulations, with each simulation receiving its own file.
    - [Sahel_1a_consolidate_cmip6_data.m](/Sahel_1a_consolidate_cmip6_data.m): used to combine CMIP6 simulations into one file.
* Post-processing:
    - [Sahel_1b_remove_buggy_models.m](/Sahel_1b_remove_buggy_models.m): should be run after downloading new data. Removes models which produced outputs off by orders of magnitude which were identified by visual examination of the data.
* Other files:
    - [concat_script.sh](/concat_script.sh): a useful, but not currently used, bash script for combining monthly netcdf files downloaded directly from ESGF into one continuous file.

# Analysis

## Analysis Scripts

### [Sahel_2_make_means.m](/Sahel_2_make_means.m)
Performs a tiered MMM and saves the results in the data folder, with the first tier (the "Model Mean") under \*\_MM.mat, and the second (the "Institution Mean", previously the "Group Mean") under \*\_GM.mat.

### [Sahel_3_analyze_means.m](/Sahel_3_analyze_means.m)
Performs bootstrapping analysis and saves the results in a folder called Analysis, labelled by simulation, years analyzed, and the number of bootstrapping iterations.

### Sahel_3...
These other files are used to create figures used in Herman et al 2020 and followup work. The figures are saved in a folder titled "figures" and subfolders by variable, labelled according to the equivalent figure number in Herman et al 2020.

## [functions](/functions)
This folder contains functions called by the data and analysis scripts.

## Other Files

* [AA_explained.m](/AA_explained.m): creates figure comparing AA simulations of NA SST to attributed NA SST when a linear GHG trend is removed from ALL simulations.
* [calc_teleconnection.m](/calc_teleconnection.m): Uses bootstrapped MMMs to calculate covariances between NARI and Sahelian precipitation and associated uncertainties.
* [make_table_1.m](/make_table_1.m): creates supplementary tables listing simulations used in followup work.
