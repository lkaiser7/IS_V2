### invasive species models source script ###
### scripts to build and run biomod2 sdms ###
### master code to use in sdm IS_analysis ###

# clear the environment, temp files, and all other variables
rm(list = ls())

##########################################
##### SET SOURCE LOCATIONS AND PATHS #####
##########################################

# set root path to source files
# rootDir<-"D:/projects/Invasives_modeling/Invasive_SDMs/"
# rootDir<-"E:/Invasive_SDMs/"
project_dirs=c("C:/Users/lkaiser-local/Desktop/Phase1_SDMs/", "E:/invasives_SDM/", "/home/pierc/projects/invasives_SDM/")
rootDir=project_dirs[min(which(dir.exists(project_dirs)))]

# set working directory to main analysis folder
setwd(rootDir)

run_type="local_HI" #local_HI global_notHI nested_HI
# select name for project and create directory
project_run<-paste0(run_type, "_models")
# set path of ongoing project run for all outputs
project_path<-paste0(rootDir, project_run, "/")
# create project folder path
dir.create(project_path, showWarnings = FALSE)

# location to save any extra or more specific outputs
outDir<-paste0(project_path, "outputs/")
# create output folder in project path
dir.create(outDir, showWarnings = FALSE)

# location of scripts and code
# codeDir<-paste0("D:/projects/Invasives_modeling/Invasive_SDMs/IS_V2/")
codeDirs=c("D:/projects/Invasives_modeling/IS_V2_repo/", paste0(rootDir, "IS_V2/"), "/home/pierc/git_repos/IS_V2/") #in order of priority
codeDir=codeDirs[min(which(dir.exists(codeDirs)))]

# location of all data
dataDir<-paste0(rootDir, "data/")
# all species data
allDir<-paste0(dataDir, "merged_data/all_data/")
# Hawaii species data
hiDir<-paste0(dataDir, "merged_data/hi_data/")
# excluding Hawaii species data
nohiDir<-paste0(dataDir, "merged_data/no_hi_data/")

# location of map data and shapefiles
mapDirs<-c(paste0(dataDir, "map_data/"), "D:/data/")
mapDir<-mapDirs[min(which(dir.exists(mapDirs)))]

# location of bioclimatic variables
if (run_type=="global_notHI"){
  bioclims_dirs=c("D:/data/global_climate/wc2.1_30s_bio_simplified/", "/home/pierc/data/global_climate/wc2.1_30s_bio_simplified/", paste0(dataDir, "bioclim_vars/")) #in order of priority
  bioclims_dir<-bioclims_dirs[min(which(dir.exists(bioclims_dirs)))]
}else{
  bioclims_dirs=c("/home/pierc/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/", 
                  "D:/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/", 
                  paste0(dataDir, "bioclim_vars/")) #in order of priority
  bioclims_dir<-bioclims_dirs[min(which(dir.exists(bioclims_dirs)))]
}

# global bioclim variables V2.1 downloaded from worldclim.org (2020)
# current (2000) bioclimatic variables @ 5 arc-min (170 km2)
fitting_bios_global<-paste0(bioclims_dir, "all_baseline/current_30s/")
# current(2000) bioclimatic variables @ 10 arc-min (340 km2)
#current_proj_bios_global<-paste0(bioclims_dir, "all_baseline/current_10min/") #not used
# ### V2.1 NOT YET AVAILABLE ###
# # future (2100) bioclimatic variables @ 10 arc-min
# future_bios<-paste0(bioclims_dir, "all_future/future_10min/")

# updated HRCM bioclims_dir ***FOR HAWAII ONLY*** (2015)
# current updated bioclimatic variabels @ 125 m
fitting_bios_HIs<-c(paste0(bioclims_dir, "all_HRCM/current_250m_redone/"), "D:/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/",
                    "/home/pierc/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/")
fitting_bios_HI<-fitting_bios_HIs[min(which(dir.exists(fitting_bios_HIs)))]

# current updated bioclimatic variables @ 500 m
#changed to new recalc values (GLOBAL AND LOCAL SEEM TO HAVE DIFFERENT UNITS!)
current_proj_bios_HI<-fitting_bios_HI #paste0(bioclims_dir, "all_HRCM/current_250m_redone/") 
# future updated bioclimatic variables @ 500 m 
future_proj_bios_HIs<-c(paste0(bioclims_dir, "all_HRCM/future_500m/"), "D:/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/")
future_proj_bios_HI<-future_proj_bios_HIs[min(which(dir.exists(future_proj_bios_HIs)))]

# GLOBAL A: allDir, fitting_bios_global, current/future_proj_bios_HI
# GLOBAL B: nohiDir, fitting_bios_global, current/future_proj_bios_HI
# LOCAL: hiDir, fitting_bios_HI, current/future_proj_bios_HI
# WEIGHTED: hiDir, fitting_bios_HI, current/future_proj_bios_HI

# select current data and bioclims to use for model approach
if (run_type=="global_notHI"){
  baseData<-nohiDir                # baseline species data (scripts 1 & 2)
  biofitRun<-fitting_bios_global     # for model fitting
  bioclim_scaling_factors=c(1/10, 1/10, 1, 1/10) #rasters were saved as integers
}else{
  baseData<-hiDir                # baseline species data (scripts 1 & 2)
  biofitRun<-fitting_bios_HI     # for model fitting
  bioclim_scaling_factors=c(1, 1, 1, 1)
}
futureData<-hiDir              # future species data (scripts 3 & 5)
biobaseRun<-current_proj_bios_HI    # for baseline projections
biofutureRun<-future_proj_bios_HI   # for future projections

##################################
##### GENERAL CONFIGURATIONS #####
##################################

# load necessary packages for this script
library("tools")
library("raster")
library("rgdal")
library("rworldmap")
library("rworldxtra")
library("maps")
library("maptools")

# set all_sp_nm = 'Clidemia_hirta' for testing and debugging
# list all 17 species names to be analyzed
all_sp_nm = c('Clidemia_hirta', 'Falcataria_moluccana', 'Hedychium_gardnerianum',
              'Lantana_camara', 'Leucaena_leucocephala', 'Melinis_minutiflora',
              'Morella_faya', 'Panicum_maximum',
              'Passiflora_tarminiana', 'Pennisetum_clandestinum', 'Pennisetum_setaceum',
              'Psidium_cattleianum', 'Setaria_palmifolia','Schinus_terebinthifolius',
              'Cyathea_cooperi', 'Miconia_calvescens', 'Ulex_europaeus')
# NOTE: Cyathea cooperi is the species synonym for Sphaeropteris cooperi
# NOTE: Passiflora tarminiana is a species synonym of Passiflora mollisima

# # Phase 1 Select Species
# all_sp_nm = c('Clidemia_hirta', 'Lantana_camara', 'Pennisetum_clandestinum', 'Psidium_cattleianum')

# # large species files to run separately
# all_sp_nm = c('Miconia_calvescens', 'Ulex_europaeus')

# # create a subset of species to run if needed (run individually for global)
# sp_sub<-c('Clidemia_hirta')
# all_sp_nm<-sp_sub

# load map data for global extent
world_map<-getMap(resolution = "high")
# load shapefile for Hawaii extent
hawaii_map<-readOGR(paste0(mapDir, "Main_Hawaiian_Islands_simple3.shp"))

# set projection to be the same for all mapping
coordSys<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'
# add coordinate reference system for projections to be same 
projection(world_map)<-coordSys
projection(hawaii_map)<-coordSys

# set map and scale (Hawaii or Global) to use for project run
if (run_type=="global_notHI"){
  map_scale<-"Global"
  map_to_use<-world_map
}else{
  map_scale<-"Hawaii"
  map_to_use<-hawaii_map
}
# list global extent from world_map
all_ext<-extent(world_map)
# create local extent for Hawaii
hi_ext<-extent(hawaii_map)

# set crop extent for project run
if (run_type=="global_notHI"){
  crop_ext<-all_ext
}else{
  crop_ext<-hi_ext
}
#############################
##### MODELLING OPTIONS #####
#############################

# select BIOMOD2 models to run (ANN, CTA, FDA,  GAM, GBM, GLM, MARS, MAXENT, RF, SRE)
models_to_run = c("MAXENT.Phillips", "GBM")  #("GAM", "GBM", "GLM", "MAXENT", "RF")
# select model evaluation methods (KAPPA, ROC, TSS)
eval_stats = c("ROC", "KAPPA", "TSS") 
# select environmental variables for models
env_var_files = c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
# create vector with bioclimatic variable names without the file extension (.tif)
var_names<-unlist(file_path_sans_ext(env_var_files))

# choose whether to plot graphs (T) or not (F)
plot_graphs = TRUE
# plotting options depending on if server or not
useRasterDef = TRUE
interpolateDef = FALSE
# apply fixes for large (T) or small (F) models to solve memory issues (script 3b)
apply_biomod2_fixes = T
# choose whether to overwrite past results (T) or not (F)
overwrite = T
# select number of computer cores for processing (max = 32)
cpucores = 5

if (cpucores==1){
  parallel_run = F
}else{
  parallel_run = T
}
### MAIN SCRIPTS ###

# RUN 'T' FOR BOTH BASELINE AND FUTURE
# script 1: to run model fitting (T) or not (F)
EM_fitting = TRUE
# script 2: to run ensemble modeling (T) or not (F)
EM_ensemble = TRUE
# script 3: to project model results (T) or not (F)
EM_project = TRUE

### AUXILIARY SCRIPTS ###

# RUN 'T' FOR BASELINE, 'F' FOR FUTURE
# script 1a: to get variable importance and evaluation score (T) or not (F)
merge_var_imp_and_mod_eval = T
# script 1b: to graph evaluation scores/variable importance (T) or not(F)
model_fit_graphs = T
# script 2a: to graphy variable importance (T) or not (F)
create_response_curves = T

# RUN 'T' FOR BOTH BASELINE AND FUTURE
# script 4: to create raster files (T) or not (F)
raster_output_creation = T 

# RUN 'F' FOR BASELINE, 'T' FOR FUTURE
# script 5: to map analog climates (T) or not (F)
create_analog_climate_map = F
# script 6: to calculate distribution shifts (T) or not (F)
calculate_distribution_shifts = F
# script 7: to create ensemble maps (T) or not (F)
species_ensemble_maps = F

##########################################
##### SPECIFIC SCRIPT CONFIGURATIONS #####
##########################################

### EM_fitting (script 1)
# number of ensemble modeling evaluation runs (set to 10 for full runs)
NbRunEval = 10
# if the models should use response points weights or not
if (run_type=="nested_HI"){
  useYweights = T
}else{
  useYweights = F 
}
# consider PAs outside of a species climate envelope (T) or not (F)
PseudoAbs_outside_CE = FALSE
# set PA density that is equal to point density within surveyed areas
dens_PAs_outside_CE = 1 
# select number of repetitions for PA selections
PA.nb.rep = 10
# select number of PAs to determine point density
number_of_PAs = 1 #Using 1 to 1, based on recommendation from https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/j.2041-210X.2011.00172.x
#if less than 100, will use value to determine total PA as number_of_PAs * n of presences, if larger, will apply actual number 
# candidate points to use only if PAs_outside_CE = F, if == 0, will use number_of_PAs   
candidatePAperPA = PA.nb.rep*5  # overridden if PseudoAbs_outside_CE = T
# strategy for selecting PAs (disk, random, sre, user.defined)
PA.strategy = "random" #using random, as recommended by https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/j.2041-210X.2011.00172.x
# set 100m equivalence distance from actual data points
equiv_100m = 0.0009430131
# set 20 km minimum distance from actual data points
if (run_type=="global_notHI"){
  PA.dist.min = 4*equiv_100m 
  thinning_dist_km=4 #also set presence thinning distance
}else{
  PA.dist.min = 1*equiv_100m 
  thinning_dist_km=1 #also set presence thinning distance
}

# to run the full models (T) or not (F)
do.full.models = TRUE

### EM_ensemble (script 2)
# sets the minimum scores to exclude models when building ensembles
eval.metric.threshold = rep(0.5, length(eval_stats)) 

### EM_project (script 3)
# select baseline (1) or future (4) projections
baseline_or_future = 1
# choose to save clamping mask (T) or not (F)
clampingMask = FALSE
# to keep clamping Mask = T saved to hard drive (T) or not (F)
memory = TRUE 
# assign projected climate data set for baseline scenario
if (baseline_or_future == 1) {
  clim_data = biobaseRun
  proj_nm = 'baseline'}
# assign projected climate data set for future scenario
if (baseline_or_future == 4) {
  clim_data = biofutureRun 
  proj_nm = 'future'}
# temporary folder for files during processing to avoid memory errors
dir_for_temp_files<-paste0(rootDir, project_run, "/temp/", baseline_or_future, "/") 
# conditions if applying fixes to BIOMOD2 (script 3b)
if (apply_biomod2_fixes) {
  # name model run based on scenario 
  maxentWDtmp = paste0("maxentWDtmp_", baseline_or_future)
  # create temporary directory file
  dir.create(dir_for_temp_files, showWarnings = F, recursive = T)
}

### raster_output_creation (script 4)
# number of projections to create raster - 1 for baseline, 2 for both
projections_to_run = 1 
# type of ensemble configurations for multiple species maps
spp_ensemble_type = "wmean" 
# for raster creation and shifted calculations for mapping
spp_ensemble_eval_stats = c('ROC', 'TSS', 'KAPPA')
# for raster creation and shifted calculations
comp_projects = c('baseline', 'future')
# for raster creation 
plot_spp_ensemble_CV = TRUE 
# for raster creation for future runs only **depreciated**
masked_spp_ensemble_map = FALSE

### create_analog_climate_map (script 5)
#create_analog_climate_map configurations for baseline (1) or future(4)
toCompareWithCurrentClimate = 1 

### calculate_distribution_shifts (script 6)
# km resolution for shifted calculations
model_resolution = 0.5

###########################
##### RUNNING SCRIPTS #####
###########################

# create log for data input for initial run
if (baseline_or_future == 1) {
  # get processor ID for R session 
  worker = paste0(Sys.Date(), "_worker", Sys.getpid())
  # create a text file name for specific ongoing processing session
  txt_nm = paste0(project_path, "data_input_log_", worker, ".txt")
  # write log file in species directory and sink console outputs to log
  sink(file(txt_nm, open = "wt"))
  # print sign posting for ongoing project run processing
  cat('\n', project_run, 'started on ', date(), '\n') 
  
  # list data used to keep record of inputs per run
  cat('\n Inputs used:', '\n')
  cat('species:'); print(all_sp_nm)
  cat('species baseline data:', baseData, '\n')
  cat('species future data:', futureData, '\n')
  cat('bioclimatic variables:', var_names, '\n')
  cat('fitting bioclim data:', biofitRun, '\n')
  cat('baseline bioclim data:', biobaseRun, '\n')
  cat('future bioclim data:', biofutureRun, '\n')
  cat('map and crop extent:', map_scale, '\n')
  cat('selected models:', models_to_run, '\n')
  cat('selected evaluation statistics:', eval_stats, '\n')
  cat('number of evalutations and repetitions:', NbRunEval, '&', PA.nb.rep, '\n')
  cat('response points (Y)weights used?', useYweights, '\n')
  
  # add any additional notes for project run if needed
  cat('\n', 'additional notes: ')
  
  # reset sink from text file to console output
  sink(NULL)
}

# start the clock to calculate processing time
ptmStart<-proc.time()

# run model fitting, ensemble models, and projections based on settings above
if (EM_fitting){  # 1 - run fitting code
  source(paste0(codeDir,"1_model_fitting.R")) 
}
if (EM_ensemble){  # 2 - run ensemble code
  source(paste0(codeDir,"2_mod_ensembles.R")) 
}
if (EM_project){  # 3 - run projection code
  source(paste0(codeDir,"3_em_projection.R"))   
  #source(paste0(codeDir,"3d_convert_grd_to_tiff.R")) #convert huge grd files to compressed tiff   
}

############################################
# auxiliary scripts based on settings above
if (merge_var_imp_and_mod_eval) { # 1a - variable importance/evaluation statistics
  source(paste0(codeDir, "aux_output_scripts/1a_mfit_tables.R"))}
if (model_fit_graphs) { # 1b - evaluation statiscs/variable importance graphs
  source(paste0(codeDir, "aux_output_scripts/1b_mfit_graphs.R"))}
if (create_response_curves) { # 2a - response curves
  source(paste0(codeDir, "aux_output_scripts/2a_resp_curves.r"))
  source(paste0(codeDir, "aux_output_scripts/2c_mean_resp_curves.R"))
}
source(paste0(codeDir, "aux_output_scripts/expert_maps.R"))
source(paste0(codeDir, "aux_output_scripts/summary_plots.R"))
source(paste0(codeDir, "aux_output_scripts/global_vs_local_model_eval.r"))

if (raster_output_creation) { # 4 - create output rasters
  source(paste0(codeDir,"4_raster_output.R"))}

############################################
# future projection scripts

if (create_analog_climate_map) { # 5 - map analog climates
  source(paste0(codeDir,"5_analog_climate.R"))}
if (calculate_distribution_shifts) { # 6 - calcuate distribution shift
  source(paste0(codeDir,"6_range_shifts.R"))}
if (species_ensemble_maps) { # 7 - map species ensembles
  source(paste0(codeDir,"7_ensemble_maps.R"))}

# stop the clock and calculate processing time to run all scripts and code
ptmEnd = proc.time() - ptmStart 
# store elapsed time per species 
p_time = as.numeric(ptmEnd[3])/length(all_sp_nm) 
# convert seconds into minutes
p_time = p_time/60 
# report processing time and show finished processing
cat('\n','It took ', p_time, "minutes (on average) to model each species with",
    length(models_to_run), "model types") 

# delete temp files
source(paste0(codeDir,"8_delele_tmp_files.R"))

#########################
### END SOURCE SCRIPT ###
#########################

