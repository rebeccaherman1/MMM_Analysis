# MMM_Analysis
Multi-Model Mean Analysis of Climate Simulations for Detection and Attribution

# TODO List
* bring the historical pr and ts files over
* add a T field to all data files, and modify script to do such
* clean Sahel_2 and also turn it into a function for better bootstrapping.
* fix off-by-one error in CI (Fig5)

# Accessing the Data
Observed JAS area-averaged Sahel precipitation (from 
[GPCC](/model_output/historical_precipitation.mat) and from [CRU](/model_output/CRU_data.mat))
and SST indices (from [ERSST](/data/obs_Jul-Sep_ERSST.mat)) are located in the [data](/data) folder. Seasonal and 
area averaged simulated data, labelled by institution, model, and run,  are located in 
the same folder, with file names ending in "\_all.mat." There is a different file for 
each unique experiment type and data source, described below.  
* historical experiments (radiative forcings)
    - h: CMIP5 historical "ALL"  
    - a: CMIP5 historical Anthropogenic Aerosols "AA"  
    - n: CMIP5 historical Natural Forcings "NAT"  
    - g: CMIP5 historical Greenhouse Gases "GHG"  
* v: CMIP6 "Vanilla" AMIP simulations (forced with SST and no radiative forcings)
* r: AMIP+RAD (SST + ALL radiative forcings). "r" includes data from all files below if the runs start by 1901. 
    - e, ERA: ERA20CM   
         - Currently, "e" uses 1901-2003, while "ERA" has the entire length of the time series and 
    an additional variable T keeping track of that.   
    - a6: CMIP6 AMIP   
        - The files called 1901 include fewer simulations than those labelled 1950, reflecting 
    the start times of those simulations.  
    - p: NOAA PSL-FACTS
    - amip: CMIP5 AMIP (from 1950)  

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
