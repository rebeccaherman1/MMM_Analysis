# MMM_Analysis
Multi-Model Mean Analysis of Climate Simulations for Detection and Attribution

# TODO List
* bring the historical pr and ts files over
* add a T field to all data files, and modify script to do such
* run the save data script on the NOAA database files, as well as on the CMIP6 "vanilla"
  amip runs.
* combine AMIP+ALL files: CMIP6, ERA20CM, NOAA
* add additional models to the umbrella file -- make sure it has all CMIP6 and NOAA models

# Accessing the Data
Observed JAS area-averaged Sahel precipitation (from 
[GPCC](/model_output/historical_precipitation.mat) and from [CRU](/model_output/CRU_data.mat))
and SST indices are located in the [model_output](/model_output) folder. Seasonal and 
area averaged simulated data, labelled by institution, model, and run,  are located in 
the same folder, with file names ending in "\_all.mat." There is a different file for 
each unique experiment type and data source, described below.  
h: CMIP5 historical "ALL"  
a: CMIP5 historical Anthropogenic Aerosols "AA"  
n: CMIP5 historical Natural Forcings "NAT"  
g: CMIP5 historical Greenhouse Gases "GHG"  
amip: CMIP5 AMIP (from 1950)  
e, ERA: ERA20CM (SST + ALL radiative forcings)  
    - Currently, "e" uses 1901-2003, while "ERA" has the entire length of the time series and 
    an additional variable T keeping track of that.   
a6: CMIP6 AMIP (includes radiative forcings)  
    - The files called 1901 include fewer simulations than those labelled 1950, reflecting 
    the start times of those simulations.  
v: CMIP6 "Vanilla" AMIP simulations (no radiative forcings)

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
