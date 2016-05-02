### HOW TO DO IT YOURSELF ### 

# a basic guide to using the IS_V2 GitHub repository #
# steps to build and run SDMs in the BIOMOD2 package #
# directions for using the 0_source_script & scripts #

# Welcome to the GitHub repository for Invasive Species Distribution Modeling in R!
# This archive of code contains the scripts necessary to run comparisons of current 
# and future climate scenarios for mulitple species. Here are basic directions how 
# to use this repository of scripts and build your own species distribution models.

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

# 3. set the location for your code, species data, map data, and climate data directories
codeDir   # line 31
dataDir   # line 34
mapDir    # line 51
bioclims  # line 54
# this will be the location where you have saved all these files 

# 4. change any basic model configurations you see fit (lines 73-77)
baseData      # baseline species data (scripts 1 & 2)
futureData    # future species data (scripts 3 & 5)
biofitRun     # for model fitting
biobaseRun    # for baseline projections
biofutureRun  # for future projections
# here you can select certain inputs based on your desired outcome
# Other inputs may be changed as well to change projected outcomes
# EX: crop_ext, models_to_run, eval_stats, env_var_files, ETC.

# 5. decide which scripts to run (T) or not (F) per scenario (lines 155-175)
# RUN 'T' FOR BOTH BASELINE AND FUTURE
EM_fitting   # script 1
EM_ensemble  # script 2
EM_project   # script 3
# RUN 'T' FOR BASELINE AND 'F' FOR FUTURE
merge_var_imp_and_mod_eval  # script 1a
model_fit_graphs            # script 1b
create_response_curves      # script 2a
# RUN 'F' FOR BASELINE AND 'T' FOR FUTURE
raster_output_creation         # script 4
create_analog_climate_map      # script 5
calculate_distribution_shifts  # script 6
species_ensemble_maps          # script 7

# 6. select number of runs/iterations and repetitions
nbRubEval  # line 181
PA.nb.rep  # line 187
# this will be done per each model per species

# 7. select which climate scenario to run (line 129)
baseline_or_future
# select baseline (1) or future (4) projections

# 8. run the entire script once with the baseline settings (steps 5 & 7)

# 9. change to future settings and rerun the entire script (steps 5 & 7)

# By running the main configuration file (0_source_script) twice, once for baseline and 
# once for future, the main and auxiliary scripts should create two climate scenarios 
# to compare the distribution and range changes of the selected invasive species.
# Outputs will be found in your created project_run directory folder.

### END OF BASIC DIY TUTORIAL ###