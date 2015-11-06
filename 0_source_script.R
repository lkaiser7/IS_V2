### invasive species models source script ###
### scripts to build and run biomod2 sdms ###
### master code to use in sdm IS_analysis ###

# clear the environment and all other variables
rm(list = ls())

################################
##### SET SOURCE LOCATIONs #####
################################

# set root path to source files
rootDir<-"Y:/PICCC_analysis/IS_analysis/"
# set working directory to main analysis folder
setwd(rootDir)

# create necessary paths to folders from root
# location of scripts and code
codeDir<-"C:/Users/Lauren/Dropbox/GitHub/IS_V2/"

# location of all data
dataDir<-paste0(rootDir, "data/")
# global GBIF data
gbifDir<-paste0(dataDir, "gbif_data/")
# regional BISON data
bisonDir<-paste0(dataDir, "bison_data/")
# local Hawaii data
localDir<-paste0(dataDir, "local_data/")
# select current data to use for project run
dataRun<-gbifDir

# location of bioclimatic variables
bioclims<-paste0(dataDir, "bioclim_vars/")
# current (2000) bioclimatic variables @ 30 arc-sec
fitting_bios<-paste0(bioclims, "current_10km/")
# current(2000) bioclimatic variables @ 10 arc-min
current_bios<-paste0(bioclims, "current_18.5km/")
# future (2100) bioclimatic variables @ 10 arc-min
future_bios<-paste0(bioclims, "future_18.5km/")

##################################
##### GENERAL CONFIGURATIONS #####
##################################

# load necessary packages for this script
library("tools")
library("raster")
library("rworldmap")
library("rworldxtra")

# select name for project and create directory
project_run<-"global_test_run2"
# set path of ongoing project run for all outputs
project_path<-paste0(rootDir, project_run, "/")
# create project folder path
dir.create(project_path, showWarnings = FALSE)

# location to save any extra or more specific outputs
outDir<-paste0(project_path, "outputs/")
# create output folder in project path
dir.create(outDir, showWarnings = FALSE)

# list all species names to be analyzed
all_sp_nm = c('Clidemia_hirta', 'Falcataria_moluccana', 'Hedychium_gardnerianum', 
              'Lantana_camara', 'Leucaena_leucocephala', 'Melinis_minutiflora', 
              'Miconia_calvescens', 'Morella_faya', 'Panicum_maximum', 
              'Passiflora_tarminiana', 'Pennisetum_clandestinum', 'Pennisetum_setaceum', 
              'Psidium_cattleianum', 'Setaria_palmifolia','Schinus_terebinthifolius', 
              'Cyathea_cooperi', 'Ulex_europaeus')
# NOTE: Cyathea cooperi is the species synonym for Sphaeropteris cooperi

# set projection for all mapping **NEEDS WORK**
coordSys<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'
# load map data for global extent
world_map<-getMap(resolution = "high")
# list global extent from world_map
all_ext<-extent(-180, 180, -90, 90)
# create local extent for Hawaii
hi_ext<-extent(-160, -155, 15, 25)
# set crop extent for project run
crop_ext<-all_ext

#############################
##### MODELLING OPTIONS #####
#############################

# select BIOMOD2 models to run (ANN, CTA, FDA,  GAM, GBM, GLM, MARS, MAXENT, RF, SRE)
models_to_run = c("GBM")  #("GAM", "GBM", "GLM", "MAXENT", "RF")
# select model evaluation methods (KAPPA, ROC, TSS)
eval_stats = c("ROC") 
# select environmental variables for models
env_var_files = c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
# create vector with bioclimatic variable names without the file extension (.grd)
var_names<-unlist(file_path_sans_ext(env_var_files))

# choose whether to plot graphs (T) or not (F)
plot_graphs = FALSE
# apply fixes for large (T) or small (F) models to solve memory issues
apply_biomod2_fixes = TRUE
# choose whether to overwrite past results (T) or not (F)
overwrite = FALSE
# turn on (T) or off (F) parallel multi-instance computing 
paralelize = FALSE
# select number of computer cores for processing (max = 32)
cpucores = 20

# to run model fitting (T) or not(F)
EM_fit = TRUE
# to run ensemble modelling (T) or not (F)
EM_ensemble = TRUE
# to project model resutls (T) or not (F) 
EM_project = TRUE

# NA? optional scripts (options below)
merge_all_var_importance_and_model_eval = F
model_eval_graph = F
var_importance_graph = F
create_response_curves = F
create_analog_climate_map = T
raster_output_creation = T #additional/ overridding configurations within files
distribution_shift_calculations = T 
spp_ensemble_maps = T

#################################
##### SCRIPT CONFIGURATIONS #####
#################################

# number of ensemble modelling evaluation runs
NbRunEval = 4
# consider PAs outside of a species climate envelope (T) or not (F)
PseudoAbs_outside_CE = FALSE
# set PA density that is equal to point density within surveyed areas
dens_PAs_outside_CE = 1 
# select number of repetitions for PA selections
PA.nb.rep = 4
# select number of PAs to determine point density
PA.nb.absences = 10000 
# candidate points to use only if PAs_outside_CE = F, if == 0, will use PA.nb.absences   
candidatePAperPA = 200  # overridden if PseudoAbs_outside_CE = T
# strategy for selecting PAs (disk, random, sre, user.defined)
PA.strategy = "random" 
# set 100m equivalence distance from actual data points
equiv_100m = 0.0009430131
# set 20 km minimum distance from actual data points
PA.dist.min = 200*equiv_100m 
# to run the full models (T) or not (F)
do.full.models = TRUE

# sets the minimum scores to exclude models when building ensembles
eval.metric.threshold = rep(0.5, length(eval_stats)) 

# configure projections by island (T) or archipelago (F) to avoid memory issues
proj_by_island = FALSE 
# select baseline (1) or future (4) projections
baseline_or_future = 4
# choose to save clamping mask (T) or not (F)
clampingMask = FALSE
# to keep clamping Mask = T saved to hard drive (T) or not (F)
memory = TRUE 

#create_analog_climate_map configurations for baseline (1) or future(4)
toCompareWithCurrentClimate = 4

#raster_output_creation configurations for multiple species maps
spp_ensemble_type = "wmean" 
# for raster creation and shifted calculations for mapping
spp_ensemble_eval_stats = c('ROC')
# for raster creation adn shifted calculations
comp_projects = c('baseline', 'future')
# for raster creation 
plot_spp_ensemble_CV = TRUE 
# for raster creation
masked_spp_ensemble_map = FALSE 
# km resolution for shifted calculations
model_resolution = 0.5
# for shifted calculations
exclude_areas_beyond_primary_habitat = FALSE
# for multiple species maps
habitat_overlay = FALSE
# overlay of prehistorical habitat for multiple species maps (landfire BPS)
BPS = FALSE

###########################
##### RUNNING SCRIPTS #####
###########################

# plotting options depending on if server or not
useRasterDef = TRUE
interpolateDef = FALSE

# temporary folder for files during processing to avoid memory errors
dir_for_temp_files<-paste0(rootDir, 'temp/', project_run, "/", baseline_or_future, "/") 

# conditions if applying fixes to BIOMOD2
if (apply_biomod2_fixes) {
  # name model run based on scenario 
  maxentWDtmp = paste0("maxentWDtmp_", baseline_or_future)
  # create temporary directory file
  dir.create(dir_for_temp_files, showWarnings = F, recursive = T)
}

# assign projected climate data set for baseline scenario
if (baseline_or_future == 1) {
  clim_data = current_bios
  proj_nm = 'baseline'}
# assign projected climate data set for future scenario
if (baseline_or_future == 4) {
  clim_data = future_bios 
  proj_nm = 'future'}

# start the clock to calculate processing time
ptmStart <- proc.time()

# run model fitting, ensemble models, and projections based on settings above
if (EM_fit){  # runs fitting code
  source(paste0(codeDir,"1_BM2_FB_SDM_fitting_LK.r")) 
}
if (EM_ensemble){  # runs ensemble code
  source(paste0(codeDir,"2_BM2_FB_SDM_EM_creation_LK.r")) 
}
if (EM_project){  #runs projection code
  source(paste0(codeDir,"3_BM2_FB_SDM_EM_projection_byIsland_LK.r"))   
}

# auxiliary scripts based on settings above
if (merge_all_var_importance_and_model_eval) {
  source(paste0(codeDir,"1opt_merge_all_var_importance_and_model_eval.R"))}
if (model_eval_graph) {
  source(paste0(codeDir,"1opt_model_eval_graph.R"))}
if (var_importance_graph) {
  source(paste0(codeDir,"1opt_var_importance_graph.R"))}
if (create_response_curves) {
  source(paste0(codeDir,"2opt_BM2_FB_SDM_response_curves.r"))}
if (create_analog_climate_map) {
  source(paste0(codeDir,"4opt_create_analog_climate_map_LK.R"))}
if (raster_output_creation) {
  source(paste0(codeDir,"4_SDM_raster_output_creation_newBiomVer.R"))}
if (distribution_shift_calculations) {
  source(paste0(codeDir,"6_distribution_shift_calculations_newbiomod2.r"))}
if (spp_ensemble_maps) {
  source(paste0(codeDir,"7_spp_ensemble_maps.r"))}


# stop the clock and calculate processing time to run all scripts and code
ptmEnd = proc.time() - ptmStart 
# store elapsed time per species 
p_time = as.numeric(ptmEnd[3])/length(spp_nm) 
# convert seconds into minutes
p_time = p_time/60 
# report processing time and show finished processing
cat('\n','It took ', p_time, "minutes (on average) to model each species with",
    length(models_to_run), "model types") 






##### PREVIOUS WORK #####
### CHANGING SCRIPT CONFIGURATIONS ###

### SCRIPT 1 ###
# global datasets
source(paste0(codeDir, "1_gbif_data.r"))
### SCRIPT 2 ###
# local datasets
source(paste0(codeDir, "2_local_data.r"))
### SCRIPT 3 ###
# bioclimatic variables
source(paste0(codeDir, "3_bioclim_vars.r"))
### SCRIPT 4 ###
# presence and absence data
source(paste0(codeDir, "4_pres_abs.r"))
### SCRIPT 5 ###
# biomod2 modeling 
source(paste0(codeDir, "5_biomod_sdm.r"))
### SCRIPT 6 ###
# modeling outputs
source(paste0(codeDir, "6_model_outputs.r"))