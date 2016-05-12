### HOW TO DO IT YOURSELF ### 

# a basic guide to using the IS_V2 GitHub repository #
# steps to build and run SDMs in the BIOMOD2 package #
# directions for using the 0_source_script & scripts #

# Welcome to the GitHub repository for Invasive Species Distribution Modeling in R!
# This archive of code contains the scripts necessary to run comparisons of current 
# and future climate scenarios for mulitple species. Here are basic directions how 
# to use this repository of scripts and build your own species distribution models.

# NOTE: THIS SCRIPT WILL NOT RUN! 
# THESE ARE ONLY INSTRUCTIONS ON HOW TO USE THIS IS_V2 REPOSITORY

###################
# LIST OF SCRIPTS #
###################

### MAIN SCRIPTS ###

# 0_source_script.R
# 1_model_fitting.R
# 2_mod_ensembles.R
# 3_em_projection.R

### AUXILIARY SCRIPTS ###

# 1a_mfit_tables.R
# 1b_mfit_graphs.R
# 2a_resp_curves.R
# 3a_project_mod.R
# 4_raster_output.R
# 5_analog_climate.R
# 6_range_shift.R
# 7_ensemble_maps.R

##############
# DIRECTIONS #
##############

# The main configuration file (0_source_script) holds all the inputs for the other scripts.
# Any of these settings and configurations can be changed in that main script to reflect
# a certain desired output.  There are a few main lines that should generally be changed.
# Setting source locations of all data and scripts will vary by each individual computer.
# These basic directions use select invasive species data available to download as a demo
# to compare baseline and future projected ranges of the distributions for these species.

### CHANGES TO 0_source_script.R FILE ###

# 1. set your source location (root directory) to your individual computer (line 13)
rootDir
# this will be unique to each individual computer with the downloaded data

# 2. assign your project run a name (line 18)
project_run
# this will be unique to the individual analysis currently running
# it will also create a folder with all the modeling outputs in it

# 3. set the location for your code, species data, map data, and climate data directories
codeDir   # line 30
dataDir   # line 33
mapDir    # line 50
bioclims  # line 53
# this will be the location where you have saved all these data files 

# 4. change any basic model configurations you see fit (lines 71-76)
baseData      # baseline species data (scripts 1 & 2)
futureData    # future species data (scripts 3 & 5)
biofitRun     # for model fitting
biobaseRun    # for baseline projections
biofutureRun  # for future projections
# here you can select certain inputs based on your desired outcome
# Other inputs may be changed as well to change projected outcomes
# EX: crop_ext, models_to_run, eval_stats, env_var_files, ETC.

# 5. select which species will be used in specific project run
all_sp_nm  # line 91
# use a subset of species names to run less species (i.e. sample data)
sp_sub     # line 101

# 6. decide which scripts to run (T) or not (F) per scenario (lines 155-185)
# RUN 'T' FOR BOTH BASELINE AND FUTURE
EM_fitting   # script 1
EM_ensemble  # script 2
EM_project   # script 3
raster_output_creation      # script 4
# RUN 'T' FOR BASELINE AND 'F' FOR FUTURE
merge_var_imp_and_mod_eval  # script 1a
model_fit_graphs            # script 1b
create_response_curves      # script 2a
# RUN 'F' FOR BASELINE AND 'T' FOR FUTURE
create_analog_climate_map      # script 5
calculate_distribution_shifts  # script 6
species_ensemble_maps          # script 7

# 7. select number of runs/iterations, weights, repetitions, and pseudo-absences
nbRubEval       # line 191
useYweights     # line 193 (for nested models if using global model suitability scores)
PA.nb.rep       # line 199
PA.nb.absences  # line 200
# this will be done per each model per species

# 8. select which climate scenario to run (line 219)
baseline_or_future
# select baseline (1) or future (4) projections

# 8. run the entire '0_source_script' once with the baseline settings (steps 6 & 8)

# 9. change to future settings and rerun the entire '0_source_script' (steps 6 & 8)

# By running the main configuration file (0_source_script) twice, once for baseline and 
# once for future, the main and auxiliary scripts should create two climate scenarios 
# to compare the distribution and range changes of the selected invasive species.
# Outputs will be found in your created project_run directory folder. The main source
# script (0_source_script) has been set to run sample species and data as an example
# with the data provided saved to your own personal computer as the location (rootDir).

# for additional question please contact Lauren R. Kaiser at lkaiser7@hawaii.edu 

### END OF BASIC DIY TUTORIAL ###