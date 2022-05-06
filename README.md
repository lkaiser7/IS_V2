# IS_V2
Invasive Species Distribution Modeling in Hawaii

A series of scripts to model global, regional, and local distributions of species using the R Statistcal Programming Software with the purpose to identify areas in Hawaii most vulnerable to invasive speces (see HOWTO_DIY file for instructions)

ABOUT THE SCRIPTS:
### 0_source_script_SAMPLE ###
sample of source code that runs all other scripts with various input variables based on IS_sample_data.zip file

0_source_script_MASTER.R more recent main interface with user selections to run species distribution models

0_source_script_MASTER_globalrun.R example for selected inputs for global runs

0_source_script_MASTER_localrun.R example for selected inputs for local runs

### 1_model_fitting ###
BIOMOD2 invasive species data formatting, options, modeling, and evaluation 

### 1a_mfit_tables ###
tables to save model fitting outputs of evaluation scores and variable importance values

### 1b_mfit_graphs ###
graphs of model fitting evaluation scores and variable importance values

### 1c_mfit_mean_tables ###
variable importance by model type per species

### 2_mod_ensembles ###
ensemble model creation from model fitting outputs

### 2a_resp_curves ###
response curve graphs from ensemble models
NOTE: there is no "2b" script

### 2c_mean_resp_curves ###
response curve graphs for mean ensemble per model type and species

### 3_em_projection ###
projections of species distributions from fitted models and created ensembles

### 3a_project_mod ###
modification to run for proper ensemble model projections

### 4_raster_output ###
creation of rasters (tiff) files and images of mapped species distributions

### 5_analog_climate ###
map creation of climate envelop showing areas with same climatic conditions

### 6_range_shifts ###
calculation of species range shifts from baseline to future runs

### 7_ensemble_maps ###
mapped species distribution from ensembles showing change in scenarios

### 8_delete_temp_files ###
short code for temporary files to be deleted to free up space and memory if needed
