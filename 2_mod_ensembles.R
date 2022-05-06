### BIOMOD2: ensemble models ###
### see source script inputs ###
### step 2: 

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)
library(biomod2)
library(stringr)

# set species name to Miconia_calvescens for debugging (all_sp_nm[7])
# sp_nm<-"Miconia_calvescens"

# reset species counter to 1
sp_nm = all_sp_nm[1] 

# initialize snowfall parallel computing function
sp_parallel_run = function(sp_nm){  
  # load necessary packages
  library(biomod2)
  library(stringr)
  require(snowfall)
  
  # convert species name to character object (in case the species are numbered)
  sp_nm = as.character(sp_nm) 
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  # create new folder with species name
  dir.create(paste0(project_path, sp_dir), showWarnings = FALSE) 
  
  # get processor ID for R session 
  worker = paste0(Sys.Date(), "_worker", Sys.getpid())
  # create a text file name for specific ongoing processing session
  txt_nm = paste0(project_path, sp_dir, proj_nm, "_ensemble_model_log_", worker, ".txt")
  # write log file in species directory and sink console outputs to log
  sink(file(txt_nm, open = "wt"))
  
  # print list start date of processing
  cat('\n', 'Started on ', date(), '\n') 
  # list initial processing time 
  ptm0<-proc.time()
  
  # record time and date prior to model fitting
  cat(format(Sys.time(), "%a %b %d %Y %X"))
  # print sign posting of ongoing model fitting
  cat('\n', sp_nm, 'ENSEMBLE CREATION:') 

  # set name of file to load workspace data from model runs
  workspace_name = paste0(project_path, sp_dir, sp_nm,"_modelfitting.RData") 
  # set name of file to save all workspace data after ensemble creation
  workspace_name_out = paste0(project_path, sp_dir, sp_nm,"_ensemblefit.RData") 
  
  # if species file already exists, skip to end, otherwise begin ensemble modeling
  if (file.exists(workspace_name_out) == FALSE | overwrite == TRUE){
    
    #####################
    ##### LOAD DATA #####
    #####################
    
    # load R enviromment from model results for first species
    load(workspace_name)
    # replace any "_" with "." in the species name
    sp_nm = str_replace_all(sp_nm,"_", ".") 

    # MODELING SUMMARY: return summary of BIOMOD2 model fitting for species
    myBiomodModelOut

    # MODEL EVALUATION: get array of all model evaluations for each species
    myBiomodModelEval<-get_evaluations(myBiomodModelOut) 
    
    # VARIABLE IMPORTANCE: get array of all model variable importance values
    get_variables_importance(myBiomodModelOut)
    
    # list models computed from model fitting
    remaining_models = myBiomodModelOut@models.computed
    
    #########################
    ### ENSEMBLE MODELING ###
    #########################
    
    # change working directory to project path to save model outputs
    setwd(project_path)
    
    # combine models and make ensemble predictions built from 1_model_fitting 
    myBiomodEM<-BIOMOD_EnsembleModeling( 
      modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
      chosen.models = remaining_models,  #vector of model runs to use
      em.by = 'all',  #how models will be combined
      eval.metric = eval_stats, #evaluation metrics to build ensemble
      eval.metric.quality.threshold = eval.metric.threshold,  #threshold to exclude models
      prob.mean = TRUE,  #estimate mean probabilities 
      prob.cv = TRUE,  #estimate coefficient of variation
      prob.ci = TRUE,  #estimate confidence interval of prob.mean
      prob.ci.alpha = 0.05,  #signficance level for estimating confidence interval
      prob.median = TRUE,  #estimate median
      committee.averaging = TRUE,  #estimate committee averaging
      prob.mean.weight = TRUE,  #estimate weighted sums
      prob.mean.weight.decay = 'proportional' ) #define relative importance of weights

    # reset working directory to root 
    setwd(rootDir)
    
    # sign-posting of completed ensemble model creation
    cat('\n', sp_nm, 'ensemble complete.') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    ########################
    ### ENSEMBLE OUTPUTS ###
    ########################
    # print enesemble modeling summary
    myBiomodEM
    
    # get evaluation scores and statistics
    get_evaluations(myBiomodEM)
    # save BIOMOD2 enesmble modeling outputs from workspace
    save("myBiomodEM", "myBiomodModelOut", "remaining_models", file = workspace_name_out)
        
    # sign-posting of completed ensemble modeling
    cat('\n all ensemble modeling complete.')
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    #Stop the clock
    # calculate total processing time
    ptm1 = proc.time() - ptm0 
    # store time elapsed and convert to minutes
    p_time = as.numeric(ptm1[3])/60 
    # report processing time to log file
    cat('\n It took', p_time, 'minutes for', sp_nm, 'ensemble modeling.')    
  }else{
    # sign-posting if file for ensemble modeling has already been created 
    cat('\n', sp_nm, 'ensemble modleing already done...')
    # indicates that this species has already been run 
  }
  # reset sink to console output
  sink(NULL)
}
# END snowfall function

########################################
##### SNOWFALL FUNCTION SCRIPT RUN #####
########################################

if (is.null(cpucores)){
  # determine total number of processors available
  cpucores = as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))  
}else{
  # determine minimum number of processors available
  cpucores = min(cpucores, as.integer(Sys.getenv('NUMBER_OF_PROCESSORS')))
}
# initialize parallel computing on minimum number of cpu cores
sfInit(parallel = parallel_run, cpus = cpucores)
# export all environmentals variables to cores for cluster programming
sfExportAll() 
# time parallel run calculation of function across cores
system.time((sfLapply(all_sp_nm, fun = sp_parallel_run)))
# remove all global environmental variables from cluster cores
sfRemoveAll()
# end parallel computing on other cpu cores
sfStop()

#################################
##### END ENSEMBLE MODELING #####
#################################

#############################
#delete temp raster files
sp_nm = all_sp_nm[1]
for (sp_nm in all_sp_nm){
  sp_nm = as.character(sp_nm) 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  temp_sp_files_to_delete<-paste0(project_path, sp_dir, "delete_temp_sp_files/", "*")
  unlink(temp_sp_files_to_delete, recursive=T, force=T) #delete previous frames
  #Loc <- "mydir"
  #system(paste0("rm -r ", temp_sp_files_to_delete))
}

temp_loc_to_delete=paste0("E:/Invasive_SDMs/global_model/temp/", "*")
unlink(temp_loc_to_delete, recursive=T, force=T) #delete previous frames
