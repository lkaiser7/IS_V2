### model fitting of BIOMOD2 ###
### see source script inputs ###
### step 1: 

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)
require(tools)

# copy maxent jar and batch files to project path for processing
file.copy(c(paste0(dataDir, "maxent/maxent.jar"), paste0(dataDir, "maxent/maxent.bat")),
          project_path, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# set sp_nm = 'Clidemia_hirta' for testing and debugging (all_sp_nm[1])

# initialize snowfall parallel computing function
sp_parallel_run = function(sp_nm) {
  # load necessary packages
  library(biomod2)
  library(raster)
  library(rgeos)
  library(randomForest)
  library(dismo)
  library(mda)
  library(stringr)
  library(tools)
  require(snowfall)
  
  # convert species name to character object (in case the species are numbered)
  sp_nm = as.character(sp_nm) 
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  # create new folder with species name
  dir.create(paste0(project_path, sp_dir), showWarnings = FALSE) 
  
  # create temporary directory per species to be deleted
  temp_sp_files_to_delete<-paste0(project_path, sp_dir, "delete_temp_sp_files/")
  dir.create(temp_sp_files_to_delete, showWarnings = FALSE)
  # set temporary directory to created temp file
  rasterOptions(tmpdir = temp_sp_files_to_delete)
  
  # print posting of temporary file location
  cat('\n temporary files to be deleted saved here:', temp_sp_files_to_delete, '\n')
 
  # get processor ID for R session 
  worker = paste0(Sys.Date(), "_worker", Sys.getpid())
  # create a text file name for specific ongoing processing session
  txt_nm = paste0(project_path, sp_dir, proj_nm, "_fitting_log_", worker, ".txt")
  # write log file in species directory and sink console outputs to log
  sink(file(txt_nm, open = "wt"))

  # print list start date of processing
  cat('\n Started on ', date(), '\n') 
  # list initial processing time 
  ptm0<-proc.time()
  
  # record time and date prior to model fitting
  cat(format(Sys.time(), "%a %b %d %Y %X"))
  # print sign posting of ongoing model fitting
  cat('\n', sp_nm, 'MODEL FITTING:') 
  
  # set name of file to save all workspace data after model run
  workspace_name = paste0(project_path, sp_dir, sp_nm, "_modelfitting.RData")
    
  # create file name for variable importance from model runs
  FileName00<-paste0(project_path, sp_dir, sp_nm, "_VariImp.csv") 
  # if species file already exists, skip to end, otherwise begin model fitting
  if (file.exists(FileName00) == FALSE | overwrite == TRUE) {     

    #####################
    ##### LOAD DATA #####
    #####################
    
    # print posting of loading spatial raster data
    cat('\n loading rasters...')
    
    # establish bioclim variable as mask as predictor variable **USE FITTING BIOS***
    predictors = raster(paste0(biofitRun, env_var_files[1])) 
    # crop predictor to selected extent 
    predictors = crop(predictors, crop_ext, projection = coordSys) 
    
    # creates raster stack all bioclimate variables to use
    for (j in 2:length(env_var_files)){ # add rest bioclim variables to "predictors"
      temp = raster(paste0(biofitRun, env_var_files[j]))
      temp = crop(temp, crop_ext, projection = coordSys)
      predictors = addLayer(predictors, temp)
    } 
    # assign names of bioclim variables to raster stack
    names(predictors)<-var_names
    # remove temporary files
    rm(j, temp)
    
    # name output image file
    tiff_name = paste0(project_path, sp_dir, sp_nm, "_env_vars_used.tif") 
    # create blank image file in working directory
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # plot bioclimatic predictors
    plot(predictors, col = rev(terrain.colors(255)), maxpixels = 100000, 
         useRaster = useRasterDef, axes = TRUE, addfun = NULL)  
    # save image file
    dev.off() 
    # print posting of .tif image saved
    cat('\n .tif image file of environmental variables used', tiff_name, 'saved. \n')
    
    # print posting of rasters loaded and saved
    cat('\n bioclimatic variables rasters loaded and saved. (Line 110)')
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # print posting of loading species point data
    cat('\n loading species data...') 
    
    # load species occurrence (presence) data 
    mySpecies<-read.csv(paste0(baseData, sp_nm, ".csv"), header = TRUE)
    # presence data handling to retain only X and Y columns
    mySpeciesOcc<-data.frame(mySpecies$decimalLongitude, mySpecies$decimalLatitude)
    # add column for presence (1) values
    mySpeciesOcc$pa<-rep(1, length(mySpeciesOcc[,1])) 
    # rename column headers
    names(mySpeciesOcc)<-c("X", "Y", "PA")
    # check header of new presence data formatting
    head(mySpeciesOcc)
    
    # sign posting for pseudo-absence handling to define points
    cat('\n defining candidate PA points... (Line 129)')
    cat('\n begin selecting points from bioclimatic predictors:')
    
    # create raster layer based on environmental response variable cells
    mySREresp<-reclassify(subset(predictors, 1, drop = TRUE), c(-Inf, Inf, 0))
    # assign shared cells from bioclims and data to '1'
    mySREresp[cellFromXY(mySREresp, mySpeciesOcc[,1:2])]<-1 
    # calculate number of real absences (generally presence-only data) - should be 0
    act_abs = dim(mySpeciesOcc[mySpeciesOcc$PA == 0,])[1] 
    # create raster area outside of all known cells with data points
    neg_mySREresp = mySREresp == 0 
    # create matrix of potential candidate pseudo-absence points
    potential_PAs = rasterToPoints(neg_mySREresp, fun = function(x){x==1}) 
    # assign number of pseudo-absences to be selected (from source code)
    n_PA_pts = PA.nb.absences 

    # extract geographic xy data from potential candidate pseudo-absences
    potential_PAs = as.data.frame(potential_PAs[,1:2]) 
    # remove any data with missing geographic information or NAs
    true_potential_PAs = potential_PAs[complete.cases(potential_PAs),] 
    # add PA column to data frame and assign 'NA' to all rows for pseudo-absences
    true_potential_PAs=cbind(true_potential_PAs, pa=rep('NA', dim(true_potential_PAs)[1],1)) 
    # rename column headers
    names(true_potential_PAs) = c('X', 'Y', 'PA') 
    # check header of new pseudo-absences data formatting
    head(true_potential_PAs) 
    
    # merge real species occurrence data with created candidate pseudo-absence points
    mySpeciesData<-data.frame(rbind(mySpeciesOcc, true_potential_PAs))
    # print posting of completed species data loading
    cat('\n species presence and absence data loaded and saved. (Line 159)') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # sign-posting for extraction of environmental data at point locations
    cat('\n extracting env vars to points...')
    
    # create matrix with cell numbers of real data from selected bioclim variables
    bioclimData<-extract(predictors, mySpeciesData[, 1:2], cellnumbers = TRUE) 
     
    if(useYweights) {
      # location of suitability scores based on data points
      suitscores<-paste0(dataDir, "suitability_scores/")
      # load suitability score raster for species from global models 
      sp_scores<-raster(paste0(suitscores, sp_nm, "_suitability_baseline_ROC_wmean.tif"))
      # take the inverse of suitability scores to give weight to absences
      inverse_scores<-1-sp_scores
      # create matrix with extracted suitability scores to use as weights
      suitability_scores<-extract(inverse_scores, mySpeciesData[, 1:2], cellnumbers = TRUE)
            
      # create data frame with species PAs, suitability scores, and bioclim points  
      SP_bioclim_suitData<-data.frame(cbind(mySpeciesData, suitability_scores, bioclimData))
      # remove any NAs or rows with missing data
      SP_bioclim_suitData<-SP_bioclim_suitData[complete.cases(SP_bioclim_suitData), ]
      # rename column headers
      names(SP_bioclim_suitData)[4:6] = c('suit_cells', 'suit_scores', 'cells')
      # set presence ponits suitability to 1 for all occurrence points
      for(score in 1:length(SP_bioclim_suitData$suit_scores)){
        # find where species occurs (PA = 1)
        if(SP_bioclim_suitData[score, 3] == 1){
          # set corresponding suitability score to 1
          SP_bioclim_suitData[score, 5]<-1
        }
      }
      # check header of new species bioclim data formatting  
      head(SP_bioclim_suitData)
      
      # create temporary vector with unique values of PA column (1 - pres, 0 - abs, NA - PA)
      jnk = c(1, 0, NA)
      # sort data so PAs are removed before any true data
      data_sort<-SP_bioclim_suitData[order(match(SP_bioclim_suitData$PA, jnk)),] 
      # remove temporary vector
      rm(jnk)    
      # identify duplicates within the cell column
      dup_data<-duplicated(data_sort$cells) 
      # calculate number of duplicate entries
      n_dups = length(dup_data[dup_data == TRUE]) 
      
      # create new data frame without duplicate cells and drop cell columns  
      SP_ALL_data<-data_sort[!dup_data, c(-4,-6)]
    } else {
      # create data frame with bioclim data points and species pres&abs data
      Species_bioclimData<-data.frame(cbind(mySpeciesData, bioclimData)) 
      # remove any NAs or rows with missing data
      Species_bioclimData<-Species_bioclimData[complete.cases(Species_bioclimData), ]
      # check header of new species bioclim data formatting
      head(Species_bioclimData)
      
      # create temporary vector with unique values of PA column (1 - pres, 0 - abs, NA - PA)
      jnk = c(1, 0, NA)
      # sort data so PAs are removed before any true data
      data_sort<-Species_bioclimData[order(match(Species_bioclimData$PA, jnk)),] 
      # remove temporary vector
      rm(jnk)    
      # identify duplicates within the cell column
      dup_data<-duplicated(data_sort$cells) 
      # calculate number of duplicate entries
      n_dups = length(dup_data[dup_data == TRUE]) 
      
      # create new data frame without duplicate cells and drop cell column  
      SP_ALL_data<-data_sort[!dup_data, -4] 
    }
      
    # print posting of primary processing of species' location predictor values
    cat('\n extration of bioclimatic predictor values per species data done. (Line 233)')

    # check header of final bioclim data formatting
    head(SP_ALL_data)
    tail(SP_ALL_data)
    # check number of presences, absences, and pseudo-absences
    table(SP_ALL_data$PA, useNA = "ifany")
    # save output as .csv file (DO NOT RUN: takes a long time and lots of memory due to size)
    # write.csv(SP_ALL_data, file = paste0(project_path, sp_dir, sp_nm, 
    #                                      "_bioclim_suit_points.csv"))
      
    # print posting of completed building of species location data and predictor values 
    cat('\n species bioclimatic data selection and duplication removal complete. (Line 245)')
    
    # name output image file
    tiff_name = paste0(project_path, sp_dir, sp_nm, "_pres_pas.tif") 
    # create blank image file in working directory
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # plot world map background
    plot(map_to_use, main = paste0(sp_nm, " Presence & Pseudo-Absence Points"))  
    # add absence points
    points(subset(SP_ALL_data, PA == "NA"), col = "blue", pch = 20)
    # add presence points
    points(subset(SP_ALL_data, PA == 1), col = "red", pch = 19)
    # add legend to map
    legend("bottomleft", c ('presence', 'pseudo-absence'), 
           pch = 16, col = c("red", "blue"))
    # save image file
    dev.off() 
    
    # print posting of saved .tif image file of PA points 
    cat('\n .tif image file of presence/absence points', tiff_name, 'saved. (Line 265) \n')
    
    # sign-posting indicating number of duplicate entries removed
    cat('\n Of', length(dup_data), 'points,', 
        n_dups, 'were within the same raster cell and removed for', sp_nm,
        'and all bioclimatic point data was extracted, mapped, and [saved].')
    if(useYweights){
      cat('\n suitability scores were also extracted and to use for model fitting.')
    }
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    #############################################
    ##### BIOMOD DATA FORMATING AND FITTING #####
    #############################################
    
    # sign-posting defining variable to be used in BIOMOD2 package
    cat('\n biomod model configuration...') 
    
    # new data frame of XY coordinates only
    myRespXY = SP_ALL_data[, 1:2]
    # new data frame for presences (1), absences (0), or PAs (NA)
    myResp<-data.frame(SP_ALL_data[, 3])
    # new data frames in biomod data input format
    if(useYweights) {
      # re-assign all 'NA' pseudo-absences to real <0> absences for weighting 
      myResp[myResp == 'NA'] <- 0
      # new data frame with bioclimatic point data only
      myExpl<-SP_ALL_data[, 5:dim(SP_ALL_data)[2]]
      # new data frame with suitability scores only for weights 
      myYweights = SP_ALL_data[, 4]
    } else {
      # re-assign all 'NA' pseudo-absences to real <NA> for formating
      myResp[myResp == 'NA'] <- NA
      # new data frame with bioclimatic point data only
      myExpl<-SP_ALL_data[, 4:dim(SP_ALL_data)[2]]
      # set response weights to NULL if not considering weighted model 
      myYweights = NULL 
    }
    
    # load BIOMOD2 data for formatting
    myBiomodData<-BIOMOD_FormatingData(
      resp.name = sp_nm,  #species name
      resp.var = myResp,  #pres/abs/pa points
      expl.var = myExpl,  #bioclim values
      resp.xy = myRespXY,  #xy coordinates
      PA.nb.rep = PA.nb.rep,  #number of PAs selections
      PA.nb.absences = n_PA_pts,  #number of PAs to select
      PA.strategy = PA.strategy,  #how to select PAs
      PA.dist.min = PA.dist.min)  #minimum distance to presences

    # sign-posting of completed biomod data formating
    cat('\n biomod data formatting complete. (Line 313)') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # set different options for selected modeling techniques from source script
    myBiomodOption<-BIOMOD_ModelingOptions(
      # GBM: Generalized Boosted Regression                                            
      GBM = list(distribution = 'bernoulli', interaction.depth = 7,  shrinkage = 0.001, 
                 bag.fraction = 0.5, train.fraction = 1, n.trees = 100, cv.folds = 10),
      # MARS: Multivariate Adaptive Regression Splines
      MARS = list(degree = 2, penalty = 2,thresh = 0.001, prune = TRUE),
      # RF: Random Forest Classification and Regression
      RF = list(do.classif = TRUE, ntree = 100, mtry = 'default', 
                max.nodes = 10, corr.bias = TRUE), 
      # MAXENT: Maximum Entropy
      MAXENT = list(maximumiterations = 100, visible = FALSE, linear = TRUE, 
                    quadratic = TRUE, product = TRUE, threshold = TRUE, hinge = TRUE, 
                    lq2lqptthreshold = 80, l2lqthreshold = 10, hingethreshold = 15, 
                    beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
                    beta_hinge = -1, defaultprevalence = 0.5))
    
    # change working directory to project path to save model outputs
    setwd(project_path)
      
    # sign-posting for ensemble model fitting
    cat('\n model fitting...')
    
    # create ensemble model from formatting data and modeling options
    myBiomodModelOut<-BIOMOD_Modeling(
      myBiomodData,  #formatted biomod data
      models = models_to_run,  #select models to run
      models.options = myBiomodOption,  #biomod options object
      NbRunEval = NbRunEval,  #number of evaluation runs*** 10
      DataSplit = 80,  #amount of data to use for training
      Yweights = myYweights,  #response points weights
      VarImport = 4,  #permuations to estimate variable importance*** 10
      do.full.models = do.full.models,  #calibrate and evaluate to all data
      models.eval.meth = eval_stats,  #evaluation metrics
      SaveObj = TRUE,  #save output object
      rescal.all.models = TRUE)  #scale to binomial GLM
    
    # reset working directory to root 
    setwd(rootDir)
    
    # return summary output of biomod model runs
    myBiomodModelOut
        
    # sign-posting of completed biomod modeling 
    cat('\n biomod data modeling complete. (Line 360)') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # return output model evaluation metrics results
    myBiomodModelEval<-get_evaluations(myBiomodModelOut) 
    # return names of model evaluations
    dimnames(myBiomodModelEval) 
    # review model evaluations
    myBiomodModelEval[eval_stats, "Testing.data",,,]
    
    # Validate selectec metrics for all tests (TSS, ROC, or KAPPA)
    if ("TSS" %in% eval_stats){
      # return variable importance for each model
      myBiomodModelEval["TSS", "Testing.data",,,] 
      # create data frame with variable importance values
      Spp_TSS<-data.frame(myBiomodModelEval["TSS", "Testing.data",,,]) 
      # assign file path name for results 
      FileName<-paste0(project_path, sp_dir, sp_nm, "_TSS.csv") 
      # create .csv file and save TSS outputs
      write.table(Spp_TSS, file = FileName, sep = ",", col.names = NA) 
    }
    if ("ROC" %in% eval_stats){
      # return variable importance for each model
      myBiomodModelEval["ROC", "Testing.data",,,]
      # create data frame with variable importance values
      Spp_ROC<-data.frame(myBiomodModelEval["ROC", "Testing.data",,,])
      # assign file path name for results 
      FileName<-paste0(project_path, sp_dir, sp_nm, "_ROC.csv")
      # create .csv file and save ROC outputs
      write.table(Spp_ROC, file = FileName, sep = ",", col.names = NA)
    }
    if ("KAPPA" %in% eval_stats){
      # return variable importance for each model
      myBiomodModelEval["KAPPA", "Testing.data",,,]
      # create data frame with variable importance values
      Spp_KAP<-data.frame(myBiomodModelEval["KAPPA", "Testing.data",,,])
      # assign file path name for results 
      FileName<-paste0(project_path, sp_dir, sp_nm, "_KAPPA.csv")
      # save .csv file and save KAPPA outputs
      write.table(Spp_KAP, file = FileName, sep = ",", col.names = NA)
    }
    
    # get variable importance of selected bioclim variables
    get_variables_importance(myBiomodModelOut) 
    # create data frame with model variable importances
    Spp_VariImp <- data.frame(get_variables_importance(myBiomodModelOut)) 
    # save .csv file of variable importances
    write.table(Spp_VariImp, file = FileName00, sep = ",", col.names = NA)
    
    # save BIOMOD modeling output from workspace
    save("myBiomodModelOut", file = workspace_name)   
  
    # sign-posting of completed biomod fitting for species 
    cat('\n all model fitting, modeling, and evaluation complete. (Line 415')
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))

    # calculate total processing time
    ptm1 = proc.time() - ptm0 
    # store time elapsed and convert to minutes
    p_time = as.numeric(ptm1[3])/60 
    # report processing time to log file
    cat('\n It took ', p_time, "minutes to model", sp_nm)
  }else{
    # sign-posting if file for variable importance has already been created 
    cat('\n fitting for', sp_nm, 'already done...') 
    # indicates that this species has already been run 
  }    
  # delete select temporary files per species once processing is finished
  unlink(temp_sp_files_to_delete, recursive = TRUE)
  
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
sfInit(parallel = TRUE, cpus = cpucores)
# export all environmentals variables to cores for cluster programming
sfExportAll() 
# time parallel run calculation of function across cores
system.time((sfLapply(all_sp_nm, fun = sp_parallel_run)))
# remove all global environmental variables from cluster cores
sfRemoveAll()
# end parallel computing on other cpu cores
sfStop()

#############################
##### END MODEL FITTING #####
#############################