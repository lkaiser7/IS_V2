### ensemble modeling curves ###
### use after 2_mod_ensmeble ### 

##################
##### SET UP #####
##################

# load necessary packages
library(biomod2)
library(stringr)
library(ggplot2)
library(reshape2)
library(plyr)

# select first models run by last name
last_model = models_to_run[length(models_to_run)]

all_dev_off=function() {
  if (length(dev.list())>0){
    for (i in dev.list()[1]:dev.list()[length(dev.list())]) {dev.off()}
  }else{
    cat("no open devices \n")
  }
}

# loop through all species
sp_nm = all_sp_nm[1]
for (sp_nm in all_sp_nm){
  # convert species names to characters
  sp_nm = as.character(sp_nm)  
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm, "_", "."), "/")
  current_sp_nm=replace_spp_names(sp_nm)
  cat('\n plotting response curves for', sp_nm)
  
  # set species name
  sp_nm0 = sp_nm
  
  # set name of file to load workspace data from model runs
  workspace_name = paste0(project_path, sp_dir, sp_nm,"_modelfitting.RData") 
  # set name of file of workspace data from ensemble creation
  if (file.exists(workspace_name)){
    # load R enviromment from model results   
    load(workspace_name)    
  }else{
    cat('\n workspace for', sp_nm, 'does not exist.')
  }
  
  ########################
  ### INITIALIZE PLOTS ###
  ########################
  # change working directory to project path to save model outputs
  setwd(project_path)
  
  # create folder in output directory for response curve plots
  rc_fold<-paste0(outDir, "response_curves/")
  dir.create(rc_fold, showWarnings = FALSE) 
  
  # select first model from runs
  model = models_to_run[1]
  # loop through all models run
  for (model in models_to_run){
    # name output image file
    tiff_name=paste0(rc_fold, sp_nm, "_response_curve_", model, "_all_vars.tif") 
    # # create blank image file
    # tiff(tiff_name, res = 300, units = "in", pointsize = 12,
    #      width = 10, height = 10, compression = "lzw")
    # if file already exists, skip to end, otherwise begin plotting
    if (file.exists(tiff_name) == FALSE | overwrite == TRUE){
      
      # load BIOMOD2 modeling outputs
      loaded_models<-BIOMOD_LoadModels(myBiomodModelOut, models = model)
      # remove???
      #loaded_models = loaded_models[1:length(loaded_models) - 1]
      loaded_models=loaded_models[grep(loaded_models, pattern = "Full")]
      
      # plot 2D response plots
      tiff_name=paste0(rc_fold, sp_nm, "_response_curve_", model, ".tif")
      
      tiff(tiff_name, res = 300, units = "in", pointsize = 12,
           width = 10, height = 10, compression = "lzw")
      response.plot2(models  = loaded_models,
                                   Data = get_formal_data(myBiomodModelOut, 'expl.var'), 
                                   show.variables = get_formal_data(myBiomodModelOut,
                                                                    'expl.var.names'),
                                   do.bivariate = FALSE,
                                   fixed.var.metric = 'mean',
                                   save.file = "no", 
                                   name = "response_curve", 
                                   ImageSize = 800, 
                                   plot = T)
      all_dev_off()

      ###################
      ### COMBO PLOTS ###
      ###################
      
      # name output image file
      tiff_name=paste0(rc_fold, sp_nm, "_response_curve_", model, "_all_vars.tif")
      # create blank image file for each bioclimatic variable
      tiff(tiff_name, res = 300, units = "in", pointsize = 12,
           width = 10, height = 10, compression = "lzw")
      response.plot2(models  = loaded_models[1],
                     Data = get_formal_data(myBiomodModelOut, 'expl.var'), 
                     show.variables = get_formal_data(myBiomodModelOut,
                                                      'expl.var.names'),
                     do.bivariate = T,
                     fixed.var.metric = 'mean',
                     save.file = "no", 
                     name = "response_curve", 
                     ImageSize = 480, 
                     plot = T)
      all_dev_off()
      
      # sign-posting of completed response curve plotting
      cat('\n response curves for ',sp_nm, model, 'finished.')  
    }else{
      # sign-posting if file for response curves has already been created 
      cat('\n response curves for ',sp_nm, model, 'already done.')  
      # all_dev_off()
    }      
  }
}

# reset working directory to root 
setwd(rootDir)

####################################
##### END RESPONSE CURVE PLOTS #####
####################################