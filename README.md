# MMM_Analysis
Multi-Model Mean Analysis of Climate Simulations for Detection and Attribution

# Accessing the Data
All observations and model output can be found in the [data](/data) folder, sorted by variable. All variables are averaged seasonally over JAS and area-averaged over the Sahel or some other ocean basin. 
Observation files are titled "observations.mat" while simulation files end with "\_all.mat". 
Observations are categorized by source, while simulations are categorized by institution, model, and run, and titled by the type of simulation (listed below).
* historical experiments (radiative forcings)
    - h: CMIP5 historical "ALL"  
    - a: CMIP5 historical Anthropogenic Aerosols "AA"  
    - n: CMIP5 historical Natural Forcings "NAT"  
    - g: CMIP5 historical Greenhouse Gases "GHG"  
    - cmip6_: same experiments, but in CMIP6.
* amip-piF: CMIP6 "Vanilla" AMIP simulations (forced with SST and preindustrial radiative forcings)
* amip-hist: CMIP6 AMIP+RAD which begin by 1850 (SST + ALL radiative forcings).  

Below are the scripts used to generate these data files. 

## [concat_script.sh](/concat_script.sh)
A bash script for concatenating montly netcdf files into one continuous file.

## [umbrella.mat](/umbrella.mat)
A handmade matfile containing model name abbreviations (MODELS) and their corresponding 
"umbrella" research institution abbreviations (ABBREV) and full names (INSTITUTION) for
all models used in Herman et al 2020 and followup work.

## Model Lists
Files ending in "\_models.txt" contain lists of the files iterated over when downloading data 
using Sahel_1_save_data...

# Analysis

## [Sahel_2_make_means.m](/Sahel_2_make_means.m)
Performs a tiered MMM and saves the first tier under \*\_MM.mat and the second under \*\_GM.mat.

## [Sahel_2a_make_plots.m](/Sahel_2a_make_plots.m)
Make the first half of Fig1 using the files created in Sahel_2.
