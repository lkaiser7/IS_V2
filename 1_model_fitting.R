### model fitting of BIOMOD2 ###
### see source script inputs ###
### step 1: 

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)
require(tools)

sp_nm=all_sp_nm[1]
# initialize snowfall parallel computing function
sp_parallel_run = function(sp_nm) {
  # load necessary packages
  library(biomod2)
  #library(raster)
  library(rgeos)
  library(randomForest)
  library(dismo)
  library(mda)
  library(stringr)
  library(tools)
  require(snowfall)
  library(terra)
  
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
  unlink(paste0(temp_sp_files_to_delete, "*"), recursive=T, force=T) #delete previous temp files if any
  file.remove(list.files(tempdir(), full.names = T, pattern = "^file")) #https://stackoverflow.com/questions/45894133/deleting-tmp-files
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
    # predictors = raster(paste0(biofitRun, env_var_files[1])) 
    # predictors=predictors*bioclim_scaling_factors[1]
    # # crop predictor to selected extent 
    # #predictors = crop(predictors, crop_ext, projection = coordSys) #LF disabled as this will take a long time
    # 
    # # creates raster stack all bioclimate variables to use
    # for (j in 2:length(env_var_files)){ # add rest bioclim variables to "predictors"
    #   temp = raster(paste0(biofitRun, env_var_files[j]))
    #   temp=temp*bioclim_scaling_factors[j]
    #   #temp = crop(temp, crop_ext, projection = coordSys) #LF disabled as this will take a long time
    #   predictors = addLayer(predictors, temp)
    # } 
    # rm(j, temp) # remove temporary files
    
    predictors=rast(paste0(biofitRun, env_var_files)) #try
    predictors=predictors*bioclim_scaling_factors #try
    # assign names of bioclim variables to raster stack
    names(predictors)<-var_names
    
    if (aggregate_rasters_for_local_runs) predictors=terra::aggregate(predictors, fact=4, fun="mean")
    #plot(predictors)
    # cat("pred correlations \n")
    # terra::layerCor(predictors, 'pearson', na.rm=T, maxcell=100000)
    
    # layerStats(predictors, 'pearson', na.rm=T)
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
    cat("will extract preds now \n")
    #cell_numbers_df<-extract(predictors[[1]], mySpeciesOcc[, 1:2], cellnumbers = TRUE) 
    cell_numbers_df<-terra::extract(predictors[[1]], mySpeciesOcc[, 1:2], cells = TRUE) 
    #View(cell_numbers_df)
    #remove repeats
    unique_cells=!duplicated(cell_numbers_df[,3])
    #n_pres_wDups=nrow(mySpeciesOcc)
    mySpeciesOcc=mySpeciesOcc[unique_cells,]
    cat("after removing raster cell duplicates, sp records went from ", length(unique_cells), " points to ", sum(unique_cells), " points \n")
    
    #memory management
    #sort( sapply(ls(),function(x){object.size(get(x))})) 
    rm(cell_numbers_df, mySpecies)
    
    #######################
    #now thin points!
    #View(mySpeciesOcc)
    mySpeciesOcc_thin=mySpeciesOcc[mySpeciesOcc$PA==1,]
    mySpeciesOcc_thin$sp=sp_nm
    library(spThin)
    if (grepl("global", project_run)){
      thinning_dist_km=4
    }else{
      thinning_dist_km=1
    }
    cat("using ", thinning_dist_km, " km thinning distance \n")
    
    thinned_dataset_full <-
      thin( loc.data = mySpeciesOcc_thin, 
            lat.col = "Y", long.col = "X", 
            spec.col = "sp", 
            thin.par = thinning_dist_km, reps = 1, 
            locs.thinned.list.return = TRUE, 
            write.files = FALSE, 
            write.log.file = FALSE)
    #View(thinned_dataset_full)
    thinned_dataset_full=thinned_dataset_full[[1]]
    thinned_dataset_full$PA=1
    names(thinned_dataset_full)=c("X","Y", "PA")
    cat("thinned presences from ", nrow(mySpeciesOcc), " points to ", nrow(thinned_dataset_full), " points \n")
    mySpeciesOcc=rbind(thinned_dataset_full, mySpeciesOcc[mySpeciesOcc$PA!=1,]) 
    
    # store number of presence points per species 
    n_Occ_pts<-length(mySpeciesOcc$PA)
    
    #memory management
    #sort( sapply(ls(),function(x){object.size(get(x))})) 
    rm(mySpeciesOcc_thin, thinned_dataset_full)
    
    
    # sign posting for pseudo-absence handling to define points
    cat('\n defining candidate PA points... (Line 131)')
    cat('\n begin selecting points from bioclimatic predictors:')
    
    # create raster layer based on environmental response variable cells
    #mySREresp<-reclassify(predictors[[1]], c(-Inf, Inf, 0))
    mySREresp<-classify(predictors[[1]], matrix(c(-Inf, Inf, 0), ncol = 3, byrow = T))
    # assign shared cells from bioclims and data to '1'
    mySREresp[cellFromXY(mySREresp, mySpeciesOcc[,1:2])]<-NA 
    # calculate number of real absences (generally presence-only data) - should be 0
    act_abs = dim(mySpeciesOcc[mySpeciesOcc$PA == 0,])[1] 
    if (number_of_PAs<100){
      n_PA_pts<-n_Occ_pts*number_of_PAs 
    }else{
      n_PA_pts<-number_of_PAs 
    }
    
    #potentially consider more robust PA selection:
    #https://www.sciencedirect.com/science/article/abs/pii/S030438001500215X
    #https://rdrr.io/cran/mopa/man/pseudoAbsences.html
    true_potential_PAs=spatSample(mySREresp, candidatePAperPA*n_PA_pts, method="random", xy=T, replace=FALSE, na.rm=T)
    #true_potential_PAs=cbind(true_potential_PAs, pa=rep('NA', dim(true_potential_PAs)[1],1)) 
    true_potential_PAs=true_potential_PAs[,-3]
    true_potential_PAs$PA=NA 
    # rename column headers
    names(true_potential_PAs) = c('X', 'Y', 'PA') 
    # check header of new pseudo-absences data formatting
    head(true_potential_PAs) 
    #View(true_potential_PAs)
    

    # merge real species occurrence data with created candidate pseudo-absence points
    mySpeciesData<-data.frame(rbind(mySpeciesOcc, true_potential_PAs))
    rm(true_potential_PAs, mySpeciesOcc, mySREresp)
    # print posting of completed species data loading
    cat('\n species presence and absence data loaded and saved. (Line 167)') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # sign-posting for extraction of environmental data at point locations
    cat('\n extracting env vars to points...')
    
    # create matrix with cell numbers of real data from selected bioclim variables
    bioclimData<-terra::extract(predictors, mySpeciesData[, 1:2], cells = TRUE) 
    bioclimData=bioclimData[, c("cell", var_names)]
    #bioclimData<-raster::extract(predictors, mySpeciesData[, 1:2], cellnumbers=TRUE) 
    
    #sort( sapply(ls(),function(x){object.size(get(x))})) 
    #rm(predictors)
    #View(bioclimData)
    
    if(useYweights) {
      # location of suitability scores based on data points
      suitscores<-paste0("global_notHI_models/", "output_rasters/")
      #suitscores<-paste0(dataDir, "suitability_scores/")
      # load suitability score raster for species from global models 
      sp_scores<-raster(paste0(suitscores, sp_nm, "_suitability_baseline_TSS_wmean.tif"))
      # create matrix with extracted suitability scores to use as weights
      suitability_scores<-raster::extract(sp_scores, mySpeciesData[, 1:2], cellnumbers = TRUE)
      # create column for weights calculated from suitability scores
      suit_weights<-matrix(0, nrow = length(suitability_scores[, 2]), ncol = 1)
      # use eq(1) from Gallien et al. (2012) to calculate Yweights
      for(ss in 1:length(suitability_scores[, 2])) {
        g_weight<-1/(1+(suitability_scores[ss, 2]/(suitability_scores[ss, 2]-1))^2)
        suit_weights[ss, 1]<-g_weight
      }
      
      # create data frame with species PAs, suitability scores, and bioclim points  
      SP_bioclim_suitData<-data.frame(cbind(mySpeciesData, 
                                            suit_weights, 
                                            bioclimData))
      # remove any NAs or rows with missing data
      # SP_bioclim_suitData<-SP_bioclim_suitData[complete.cases(SP_bioclim_suitData), ]
      SP_bioclim_suitData<-SP_bioclim_suitData[complete.cases(bioclimData), ]
      # set presence ponits suitability to 1 for all occurrence points
      for(score in 1:length(SP_bioclim_suitData$suit_weights)){
        # find where species occurs (PA = 1)
        # if(SP_bioclim_suitData[score, 3] == 1){
        if(!is.na(SP_bioclim_suitData[score, 3])){
          # set corresponding suitability score to 1
          SP_bioclim_suitData[score, 4]<-1
        }
      }
      # check header of new species bioclim data formatting  
      head(SP_bioclim_suitData)
      tail(SP_bioclim_suitData)
      
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
      if (length(dup_data)>0){
        SP_ALL_data<-data_sort[!dup_data, -5]
      }else{
        SP_ALL_data<-data_sort[, -5]
      }
      #View(SP_ALL_data)
      
    } else {
      # create data frame with bioclim data points and species pres&abs data
      Species_bioclimData<-data.frame(cbind(mySpeciesData, bioclimData)) 
      # remove any NAs or rows with missing data
      #Species_bioclimData<-Species_bioclimData[complete.cases(Species_bioclimData), ]
      Species_bioclimData<-Species_bioclimData[complete.cases(bioclimData), ] #i think before the PAs were being removed in this step
      # check header of new species bioclim data formatting
      head(Species_bioclimData)
      
      # create temporary vector with unique values of PA column (1 - pres, 0 - abs, NA - PA)
      jnk = c(1, 0, NA)
      # sort data so PAs are removed before any true data
      data_sort<-Species_bioclimData[order(match(Species_bioclimData$PA, jnk)),] 
      # remove temporary vector
      rm(jnk)    
      # identify duplicates within the cell column
      dup_data<-duplicated(data_sort$cell) 
      # calculate number of duplicate entries
      n_dups = length(dup_data[dup_data == TRUE]) 
      
      # create new data frame without duplicate cells and drop cell column  
      if (length(dup_data)>0){
        SP_ALL_data<-data_sort[!dup_data, -4] 
      }else{
        SP_ALL_data<-data_sort[, -4]
      }
    }
    
    # print posting of primary processing of species' location predictor values
    cat('\n extration of bioclimatic predictor values per species data done. (Line 247)')
    
    # check header of final bioclim data formatting
    head(SP_ALL_data)
    tail(SP_ALL_data)
    dim(SP_ALL_data)
    # check number of presences, absences, and pseudo-absences
    table(SP_ALL_data$PA, useNA = "ifany")
    # save output as .csv file (takes a long time and space due to file size)
    write.csv(SP_ALL_data, file = paste0(project_path, sp_dir, sp_nm, "_bioclim_points.csv"))
    
    # print posting of completed building of species location data and predictor values 
    cat('\n species bioclimatic data selection and duplication removal complete. (Line 270)')
    
    # name output image file
    tiff_name = paste0(project_path, sp_dir, sp_nm, "_pres_pas.tif") 
    # create blank image file in working directory
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # plot world map background
    plot(map_to_use, main = paste0(sp_nm, " Presence & Pseudo-Absence Points"))  
    # add absence points
    tmp_points=subset(SP_ALL_data, is.na(PA))
    if (dim(tmp_points)[1]>n_PA_pts){
      tmp_points=tmp_points[c(1:n_PA_pts),]
    }
    points(tmp_points, col = "blue", pch = 20)
    # add presence points
    points(subset(SP_ALL_data, PA == 1), col = "red", pch = 19)
    # add legend to map
    legend("bottomleft", c ('presence', 'pseudo-absence'), 
           pch = 16, col = c("red", "blue"))
    # save image file
    dev.off() 
    rm(map_to_use)
    
    # print posting of saved .tif image file of PA points 
    cat('\n .tif image file of presence/absence points', tiff_name, 'saved. (Line 278) \n')
    
    # sign-posting indicating number of duplicate entries removed
    cat('\n Of', length(dup_data), 'points,', 
        n_dups, 'were within the same raster cell and removed for', sp_nm,
        'and all bioclimatic point data was extracted, mapped, and [saved].')
    if(useYweights){
      cat('\n suitability scores were also extracted to use for model fitting.')
    }
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    #############################################
    ##### BIOMOD DATA FORMATING AND FITTING #####
    #############################################
    
    # sign-posting defining variable to be used in BIOMOD2 package
    cat('\n biomod model configuration...') 
    
    if (is.null(seed)){
      # store number of random seed to ensure random sampling and reproducibility 
      rng_seed<-as.numeric(round(Sys.time(), 0))
    }else{
      rng_seed=seed
    }
    # print random seed number for record
    cat('\n biomod formating and fitting using random seed', rng_seed, "\n" )
    # set seed
    set.seed(rng_seed)
    
    # new data frame of XY coordinates only
    myRespXY = SP_ALL_data[, 1:2]
    # new data frame for presences (1), absences (0), or PAs (NA)
    myResp<-data.frame(SP_ALL_data[, 3])
    # re-assign all 'NA' pseudo-absences to real <NA> for formating
    myResp[myResp == 'NA'] <- NA
    # new data frames in biomod data input format
    if(useYweights) {
      # new data frame with bioclimatic point data only
      myExpl<-SP_ALL_data[, 5:dim(SP_ALL_data)[2]]
      # new data frame with suitability scores only for weights 
      myYweights = SP_ALL_data[, 4]
    } else {
      # new data frame with bioclimatic point data only
      myExpl<-SP_ALL_data[, 4:dim(SP_ALL_data)[2]]
      # set response weights to NULL if not considering weighted model 
      myYweights = NULL 
    }
    
    if (PA.strategy=='user.defined'){ #this blocking approach was not complete as there is no clear way to subsample PAs given they are coded as NAs.
      ##############################
      ##############################
      #block CV
      library(blockCV)
      tmp=cbind(myRespXY, myResp)
      names(tmp)=c("x", "y", "occ")
      #View(tmp)
      pa_data <- sf::st_as_sf(tmp, coords = c("x", "y"), crs = 4326)
      #Let’s plot the species data using tmap package:
      library(tmap)
      tm_shape(predictors[[1]]) +
        tm_raster(legend.show = FALSE, n = 30, palette = gray.colors(10)) +
        tm_shape(pa_data) +
        tm_dots(col = "occ", style = "cat", size = 0.1, colorNA="red")
      
      sb1 <- cv_spatial(x = pa_data,
                        column = "occ", # the response column (binary or multi-class)
                        k = PA.nb.rep, # number of folds
                        size = 10000, # size of the blocks in metres
                        selection = "random", # random blocks-to-fold
                        iteration = 50, # find evenly dispersed folds
                        biomod2 = TRUE) # also create folds for biomod2
      #sb1$biomod_table
      
      #PA.table #as many columns as PA.nb.rep
      #PA.strategy= 'user.defined' 
      myBiomodData<-BIOMOD_FormatingData(
        resp.name = sp_nm,  #species name (character)
        resp.var = myResp,  #pres/abs/pa points #myResp (1 col df)
        expl.var = stack(predictors), #myExpl,  #bioclim values (df)
        resp.xy = myRespXY,  #xy coordinates #as.matrix(myRespXY) (df)
        PA.nb.rep = PA.nb.rep,  #number of PAs selections (numeric)
        PA.nb.absences = n_PA_pts,  #number of PAs to select (numeric)
        PA.strategy = PA.strategy,  #how to select PAs
        PA.dist.min = PA.dist.min,
        PA.table = sb1$biomod_table)  #minimum distance to presences
      
    }else{
      ##############################
      #format  data for biomod
      myBiomodData<-BIOMOD_FormatingData(
        resp.name = sp_nm,  #species name (character)
        resp.var = myResp,  #pres/abs/pa points #myResp (1 col df)
        expl.var = stack(predictors), #myExpl,  #bioclim values (df)
        resp.xy = myRespXY,  #xy coordinates #as.matrix(myRespXY) (df)
        PA.nb.rep = PA.nb.rep,  #number of PAs selections (numeric)
        PA.nb.absences = n_PA_pts,  #number of PAs to select (numeric)
        PA.strategy = PA.strategy,  #how to select PAs
        PA.dist.min = PA.dist.min)  #minimum distance to presences
    }
    
    
    
    file_name_out = paste0(project_path, sp_dir, sp_nm, "_BiomodData.RData")
    save("myBiomodData", file = file_name_out)   
    
    # sign-posting of completed biomod data formating
    cat('\n biomod data formatting complete. (Line 335)') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # if (!only_save_biomod_input_data){
    ###################################
    ###################################
    #model tunning
    gbm.tree.complexity=7
    gbm.learning.rate=0.001
    gbm.bag.fraction=0.5
    if (optimize_model_params){
      ##############################
      cat('\n optimizing gbm parameters \n')
      #GBM optimization
      #https://rdrr.io/cran/dismo/man/gbm.step.html
      library(dismo)
      index_of_Ps=which(SP_ALL_data$PA==1)
      index_of_PAs=which(is.na(SP_ALL_data$PA))
      chosen_PAs=sample(index_of_PAs, size=n_PA_pts*4, replace = F)
      tmp_SP_ALL_data=SP_ALL_data[c(index_of_Ps,chosen_PAs),]
      tmp_SP_ALL_data[is.na(tmp_SP_ALL_data$PA),"PA"]=0    
      #View(tmp_SP_ALL_data)
      
      if (is.null(myYweights)){
        Ywgts=rep(1, nrow(tmp_SP_ALL_data))
      }else{
        Ywgts=myYweights[c(index_of_Ps,chosen_PAs)]      
      }
      optim_gbm=gbm.step(data=tmp_SP_ALL_data, gbm.x = var_names, gbm.y = "PA", site.weights=Ywgts,
                         n.trees = 20, max.trees = 2000,
                         tree.complexity = gbm.tree.complexity, learning.rate = gbm.learning.rate, 
                         bag.fraction = gbm.bag.fraction, n.folds = 10,) #higher tree complexity leads to more stable number of trees
      cat("optimum number of GBM trees is ", optim_gbm$n.trees, "\n")
      
      #########################
      #maxent optimization
      cat('\n optimizing maxent parameters \n')
      #SDMtune or ENMeval packages
      library(ENMeval)
      #https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12261
      #http://cran.nexr.com/web/packages/ENMeval/vignettes/ENMeval-vignette.html
      #https://cran.r-project.org/web/packages/ENMeval/ENMeval.pdf
      library(rJava)
      packages=as.data.frame(installed.packages())$Package
      if ("rJava" %in% packages){
        library("rJava")
        algo_type="maxent.jar"
      }else{
        algo_type="maxnet"
        library(maxnet)
        #library(ecospat)
      }
      eval3 <- ENMevaluate(occs=tmp_SP_ALL_data[tmp_SP_ALL_data$PA==1,-3], 
                           bg=tmp_SP_ALL_data[tmp_SP_ALL_data$PA==0,-3], 
                           partitions='randomkfold', algorithm=algo_type, 
                           tune.args=list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT", 'LQP'), rm = 1:4))
      #L = linear, Q = quadratic, H = hinge, P = product and T = threshold
      #fc are data transforms, rm is regularization
      #eval3@tune.settings
      best_model=eval3@results[which(eval3@results$delta.AICc==0),]
      if (nrow(best_model)>1) best_model=best_model[1,] #if multiple models just as good, pick simpler
      fc=as.character(best_model$fc)
      maxent_best_params=data.frame(L=grepl("L", fc), Q=grepl("Q", fc), H=grepl("H", fc),
                                    P=grepl("P", fc), T=grepl("T", fc), rm=as.numeric(as.character(best_model$rm)),
                                    gbm.trees=optim_gbm$n.trees)
      #View(eval3@results)
      FileName<-paste0(project_path, sp_dir, sp_nm, "_optim_maxent_and_GBM_parameters.csv") 
      write.csv(maxent_best_params, file = FileName, row.names = F) 
    }else{
      maxent_best_params=data.frame(L=T, Q=T, H=T,
                                    P=T, T=T, rm=1,
                                    gbm.trees=1, gbm.shrinkage = gbm.learning.rate, #0.001 
                                    gbm.bag.fraction = gbm.bag.fraction) #0.75
    }
    
    
    ###################################
    ###################################
    # set different options for selected modeling techniques from source script
    myBiomodOption<-BIOMOD_ModelingOptions(
      # GBM: Generalized Boosted Regression                                            
      GBM = list(distribution = "bernoulli", interaction.depth = gbm.tree.complexity,  shrinkage = maxent_best_params$gbm.shrinkage, 
                 bag.fraction = maxent_best_params$gbm.bag.fraction, train.fraction = 1, n.trees = maxent_best_params$n.trees, cv.folds = 10,
                 n.cores = 1), # to avoid parallel problems and models failing
      # # MARS: Multivariate Adaptive Regression Splines
      # MARS = list(degree = 2, penalty = 2,thresh = 0.001, prune = TRUE),
      # # RF: Random Forest Classification and Regression
      # RF = list(do.classif = TRUE, ntree = 100, mtry = 'default', 
      #           max.nodes = 10, corr.bias = TRUE), 
      # MAXENT: Maximum Entropy
      MAXENT.Phillips = list(path_to_maxent.jar = paste0(dataDir,"maxent/"),
                             maximumiterations = 100, visible = FALSE, 
                             linear = maxent_best_params$L, quadratic = maxent_best_params$Q, 
                             product = maxent_best_params$P, threshold = maxent_best_params$T, 
                             hinge = maxent_best_params$H, 
                             lq2lqptthreshold = 80, l2lqthreshold = 10, hingethreshold = 15, 
                             betamultiplier=maxent_best_params$rm,  beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
                             beta_hinge = -1, defaultprevalence = 0.5)
    )
    
    # change working directory to project path to save model outputs
    setwd(project_path)
    
    # sign-posting for ensemble model fitting
    cat('\n model fitting...')
    
    
    if (packageVersion("biomod2")=="4.0"){
      # create ensemble model from formatting data and modeling options
      myBiomodModelOut<-BIOMOD_Modeling(
        myBiomodData,  #formatted biomod data
        models = models_to_run,  #select models to run
        #models.options = myBiomodOption,  #biomod options object
        bm.options = myBiomodOption,  #biomod options object
        #NbRunEval = NbRunEval,  #number of evaluation runs*** 10
        nb.rep = NbRunEval,  #number of evaluation runs*** 10
        #DataSplit = 80,  #amount of data to use for training
        data.split.perc = 80,  #amount of data to use for training
        #Yweights = myYweights,  #response points weights
        weights = myYweights,  #response points weights
        #Prevalence = NULL, #used to build 'weighted response weights' 
        prevalence = NULL, #used to build 'weighted response weights' 
        #VarImport = 4,  #permuations to estimate variable importance*** 10
        var.import = 4,  #permuations to estimate variable importance*** 10
        do.full.models = do.full.models,  #calibrate and evaluate to all data
        #models.eval.meth = eval_stats,  #evaluation metrics
        metric.eval = eval_stats,  #evaluation metrics
        #SaveObj = TRUE,  #save output object
        save.output = TRUE,  #save output object
        #rescal.all.models = TRUE)  #scale to binomial GLM
        scale.models = F)  #scale to binomial GLM
      
    }else{
      # create ensemble model from formatting data and modeling options
      if (.Platform$OS.type == "unix"){system("export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-amd64/")}
      if (.Platform$OS.type == "unix"){system('export PATH="/usr/lib/jvm/java-1.11.0-openjdk-amd64:$PATH"')}
      myBiomodModelOut<-BIOMOD_Modeling(
        myBiomodData,  #formatted biomod data
        models = models_to_run,  #select models to run
        models.options = myBiomodOption,  #biomod options object
        NbRunEval = NbRunEval,  #number of evaluation runs*** 10
        DataSplit = 80,  #amount of data to use for training
        Yweights = myYweights,  #response points weights
        Prevalence = NULL, #used to build 'weighted response weights' 
        VarImport = 4,  #permuations to estimate variable importance*** 10
        do.full.models = do.full.models,  #calibrate and evaluate to all data
        models.eval.meth = eval_stats,  #evaluation metrics
        SaveObj = TRUE,  #save output object
        rescal.all.models = TRUE)  #scale to binomial GLM
      
    } 
    # reset working directory to root 
    setwd(rootDir)
    
    # return summary output of biomod model runs
    myBiomodModelOut
    
    # sign-posting of completed biomod modeling 
    cat('\n biomod data modeling complete. (Line 385)') 
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
    cat('\n all model fitting, modeling, and evaluation complete. (Line 440)')
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    # calculate total processing time
    ptm1 = proc.time() - ptm0 
    # store time elapsed and convert to minutes
    p_time = as.numeric(ptm1[3])/60 
    # report processing time to log file
    cat('\n It took ', p_time, "minutes to model", sp_nm)
    # }
    sink(NULL)
  }else{
    # sign-posting if file for variable importance has already been created 
    cat('\n fitting for', sp_nm, 'already done...') 
    sink(NULL)
    unlink(txt_nm) #delete previous frames
    # indicates that this species has already been run 
  }    
  
  # delete select temporary files per species once processing is finished
  unlink(paste0(temp_sp_files_to_delete, "*"), recursive=T, force=T) #delete previous frames
  gc()
  # reset sink to console output
}
# END snowfall function


########################################
##### SNOWFALL FUNCTION SCRIPT RUN #####
########################################
library(parallel)
if (is.null(cpucores)){
  # determine total number of processors available
  cpucores = detectCores(all.tests = FALSE, logical = TRUE) #as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))  
}else{
  # determine minimum number of processors available
  cpucores = min(cpucores, detectCores(all.tests = FALSE, logical = TRUE))
}
# initialize parallel computing on minimum number of cpu cores
sfInit(parallel = parallel_run, cpus = cpucores)
# export all environmentals variables to cores for cluster programming
sfExportAll() 
# time parallel run calculation of function across cores
system.time(sfLapply(all_sp_nm, fun = sp_parallel_run))
# remove all global environmental variables from cluster cores
sfRemoveAll()
# end parallel computing on other cpu cores
sfStop()

file.remove(list.files(tempdir(), full.names = T, pattern = "^file")) #https://stackoverflow.com/questions/45894133/deleting-tmp-files

tryCatch({options(warn=-1); detach(package:rworldmap)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:rworldxtra)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:maps)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:maptools)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:spThin)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:terra)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:ENMeval)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:rJava)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:maxnet)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 
tryCatch({options(warn=-1); detach(package:ecospat)}, error = function(err) {print(paste("MY_ERROR:  ",err))}) 


# END tryCatch
#############################
##### END MODEL FITTING #####
#############################
