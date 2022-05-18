### BIOMOD2: sdm projections ###
### see source script inputs ###
### step 3:

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)

temp_loc_to_delete=paste0("E:/Invasive_SDMs/global_model/temp/", "*")
unlink(temp_loc_to_delete, recursive=T, force=T) #delete previous frames

# reset sp_nm counter to first species
sp_nm = all_sp_nm[1] 

#overwrite=T

# initialize parallel function
sp_parallel_run = function(sp_nm){
  # load packaged in snowfall
  library(biomod2)
  library(stringr)
  library(colorRamps)
  library(rasterVis)
  library(tools)
  library(ncdf4)
  library(gbm) 
  library(raster)
  
  # set options for BIOMOD2 fixes in code (if assigned TRUE in source script)
  if(apply_biomod2_fixes) {
    # set raster package options
    raster::rasterOptions(tmpdir = dir_for_temp_files, timer = T, progress = "text", todisk = T) 
    # all model projection fixes to BIOMOD2 code created by Adam Vorsino
    source(paste0(codeDir,"3a_project_mod.R")) 
  }else{
    raster::rasterOptions(tmpdir = dir_for_temp_files, timer = T, progress = "text", todisk = T) 
  }
  
  # convert species name to character object (in case the species are numbered)
  sp_nm = as.character(sp_nm) 
  
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  
  #remove previous temp results
  temp_sp_files_to_delete<-paste0(project_path, sp_dir, "delete_temp_sp_files/", "*")
  unlink(temp_sp_files_to_delete, recursive=T, force=T) #delete previous frames
  
  # get processor ID for R session 
  worker = paste0(Sys.Date(), "_worker", Sys.getpid())
  # create a text file name for specific ongoing processing session
  txt_nm = paste0(project_path, sp_dir, proj_nm, "_em_projections_log_", worker, ".txt")
  # write log file in species directory and sink console outputs to log
  sink(file(txt_nm, open = "wt"))
  
  # print list start date of processing
  cat('\n', 'Started on ', date(), '\n') 
  # list initial processing time 
  ptm0<-proc.time()
  
  # record time and date prior to model projections
  cat(format(Sys.time(), "%a %b %d %Y %X"))
  # print sign posting of ongoing model projections
  cat('\n', sp_nm, 'MODEL PROJECTION:') 
  
  # print sign posting of loading previous workspace data from ensembles
  cat('\n loading EM_fitting workspace file...')
  
  # set file name of workspace data to load from ensemble fitting
  workspace_name = paste0(project_path, sp_dir, sp_nm,"_ensemblefit.RData") 
  # load workspace from previous step
  load(workspace_name)   
  # set name of file to save workspace data from ensemble projections
  workspace_name_out = paste0(sp_nm, "_", proj_nm, "_projections.RData") 
  
  # if species file already exists, skip to end, otherwise begin ensemble projections
  if(file.exists(workspace_name_out) == FALSE | overwrite == TRUE) {
    
    #####################
    ##### LOAD DATA #####
    #####################
    
    # sign-posting of loading bioclimatic variables data
    cat('\n loading predictors...')    
    # sign-posting indicating which raster is used for projections
    cat('\n using these env files for projection raster:', env_var_files, 
        '\n from dir:', clim_data)
    
    # create raster layer from bioclimatic variables
    predictors = raster(paste0(clim_data, env_var_files[1])) 
    # crop predictors to raster resolution extent
    #predictors = crop(predictors,  crop_ext) #not cropping it now
    # assign temporary variable to count number of bioclimatic variables for runs
    jnk0 = length(env_var_files)
    
    # add all other bioclimatic variables raster layers to raster stack
    for (jj in 2:jnk0) {
      # create temporary raster layer of each new bioclim variable layer to add
      temp = raster(paste0(clim_data, env_var_files[jj]))
      # crop temporary raster layer to raster resolution extent
      #temp = crop(temp,  crop_ext) #not cropping anymore
      # add new bioclimatic variable layer to raster stack
      predictors = addLayer(predictors, temp) 
    }
    
    # create vector of bioclimatic variables names without file extensions
    var_names<-unlist(file_path_sans_ext(env_var_files)) 
    # assign names of raster stack layers to bioclimatic variable names
    names(predictors)<-var_names 
    # remove temporary variables
    rm("temp", "jnk0") 
    # return summary of raster stack of predictors
    predictors 
    
    # sign-posting noting completion of loading predictors
    cat('\n done loading predictors.')  
    
    # set projection name to scenario (for baseline or future)
    projection_name = proj_nm
    
    # define extent of each different island (extra from old code)
    Kauai = c(-159.82,-159.26, 21.84, 22.25)
    Oahu = c(-158.32, -157.62,  21.22, 21.73)
    Molokai = c(-157.34, -156.69, 21.03, 21.25)
    Lanai = c(-157.08, -156.78, 20.70, 20.92)
    Maui= c(-156.8, -155.53, 20.46, 21.05)
    Hawaii = c(-156.10,-154.74, 18.87, 20.30)
    Kahoolawe = c(-156.8, -156.51, 20.46, 20.62)
    
    # create blank list of island names for species
    spIslandNames = c("")
    
    # set location to save R workspace data after all ensemble model projections
    workspace_name_out0 = paste0(project_path, sp_dir, sp_nm, 
                                 "_", proj_nm, "_all_models_proj.RData") 
    # check if file has already been created or not 
    if (file.exists(workspace_name_out0) == FALSE | overwrite == TRUE) { 
      
      # sign-posting on ongoing plotting of predictors
      cat('\n plotting predictors...')
      
      # for baseline (current) runs
      if (baseline_or_future == 1){ 
        # save name of jpeg file to be created for baseline scenarios
        tiff_name = paste0(project_path, sp_dir, sp_nm, "_",
                           projection_name, "_env_vars_used_for_projection.tif") 
        # create blank jpeg file in directory
        tiff(tiff_name, res = 300, units = "in", pointsize = 12,
             width = 10, height = 10, compression = "lzw")
        
        # plot predictors for baseline projections
        plot(predictors, col = rev(terrain.colors(255)), maxpixels = 100000, axes = TRUE,
             useRaster = useRasterDef, addfun = NULL, interpolate = interpolateDef) 
        
        # save tiff file
        dev.off()
      }
      # sign-posting to project BIOMOD2 outputs
      cat('\n run biomod_projection...')
      
      ######################
      ##### projection #####
      ######################
      
      # change working directory to project path to save model outputs
      setwd(project_path)
      
      # for baseline projections
      myBiomodProj_baseline<-BIOMOD_Projection(
        modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
        new.env = predictors,  #new explanatory variables to project onto for future
        proj.name = projection_name,  #projection directory
        selected.models = remaining_models,  #which models to use for projections
        binary.meth = eval_stats,  #evaluation method statistics 
        compress = 'xz',  #compression format of files on hard drive
        build.clamping.mask = clampingMask,  #if clamping mask should be saved or not
        keep.in.memory = memory) #if clamping mask should be saved to hard disk or not
      
      # reset working directory to root 
      setwd(rootDir)
      
      # reclaim memory no longer in use and return memory usage summary
      gc()
      # sign-posting of completed projections
      cat('\n', sp_nm, 'projection complete.')
      # sign-posting of memory size and availability
      cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn') 
      
      # save projection R workspace environement 
      save("myBiomodProj_baseline", file = workspace_name_out0)
      # sign-posting of completed BIOMOD2 projections
      cat('\n done with running biomod_projection.')
      
      # otherwise if projection workspace already exists
    } else {
      
      # load existing R workspace
      load(workspace_name_out0)
      # sign-posting indicating models are loaded from past runs
      cat('\n', sp_nm, 'projection of individual models loaded from past run.')
    }
    
    # run BIOMOD2 fixes if TRUE in source script
    if (apply_biomod2_fixes){
      # load baseline projections manually from directory
      myBiomodProjection<-LoadProjectionManually(myBiomodProj_baseline)        
    } else {
      # use previously created baseline projections
      myBiomodProjection<-myBiomodProj_baseline
    }
    
    # sign-posting to run forecasting of ensemble models
    cat('\n run ensemble forecasting...')
    
    # set location to save R workspace data and ensemble projections
    workspace_name_out1 = paste0(project_path, sp_dir, sp_nm, 
                                 "_", proj_nm, "_all_model_ensemble_proj.RData") 
    # check if ensemble projections workspace already exists or not
    if (file.exists(workspace_name_out1) == FALSE | overwrite == TRUE) { 
      
      
      ################################
      ##### ensemble_forecasting #####
      ################################
      
      # change working directory to project path to save model outputs
      setwd(project_path)
      
      # run ensemble projections for species
      myBiomodEF <- BIOMOD_EnsembleForecasting(
        projection.output = myBiomodProjection,  #BIOMOD.projection.out from projections
        total.consensus = TRUE,  #mean of all combined model projections
        EM.output = myBiomodEM,  #BIOMOD.EnsembleModeling.out from ensemble modeling
        binary.meth = eval_stats,  #evaluation method statistics 
        keep.in.memory = memory)  #if output should be saved to hard disk or not
      
      # reset working directory to root 
      setwd(rootDir)
      
      # sign-posting of completed ensemble projections
      cat('\n', sp_nm, 'ensemble projection done.') 
      
      # save ensemble projections R workspace environement 
      save("myBiomodProjection", "myBiomodEF", file = workspace_name_out1)
      
      # remove any temporary files
      removeTmpFiles(h = 1) 
    } else {
      # sign-posting if projection for species is already done
      cat('\n', sp_nm, 'projections previously done.')
    }
    ######################
    ######################
    cat('\n', 'converting grd to tiff files.')
    for (eval_stat in spp_ensemble_eval_stats){
      # convert species name to character
      sp_nm = as.character(sp_nm)  
      # store species name as character
      sp_nm0 = sp_nm
      # replace species naming convention of "_" with "."
      sp_nm = str_replace_all(sp_nm, "_", ".")
      
      # print sign posting of ongoing modeling per species 
      cat('\n',sp_nm,'grid conversion')
      
      # store raster names
      raster_names = c("EM_suitability1", "EM_suitability2")
      # store bin names
      raster_names_bin = c("EM_BIN1", "EM_BIN2")
      
      # reset i counter to 1
      i = 1
      # loop through 2
      for (i in 1:projections_to_run){
        # store individual project name
        proj_nm = comp_projects[i]
        # store individual raster name
        raster_name = raster_names[i]
        # store individual bin name
        raster_name_bin = raster_names_bin[i]
        
        #file_name=file_name1
        replace_grd_fx=function(file_name){
          if (file.exists(file_name)){
            library(terra) #using terra package to save as raster package does not keep layer names
            #https://stackoverflow.com/questions/26763013/r-write-rasterstack-and-preserve-layer-names
            # store first file name for raster .grd file
            temp_raster = rast(file_name) 
            output_nm=sub(pattern = ".grd$", replacement = ".tif", file_name)
            if (file.exists(output_nm) == FALSE | overwrite == 1){
              terra::writeRaster(temp_raster, filename = output_nm, gdal=c(compress="LZW"), overwrite=T)
            }
            unlink(file_name, recursive=T, force=T)
            unlink(sub(pattern = ".grd$", replacement = ".gri", file_name), recursive=T, force=T)
          }
        }
        
        file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                            proj_nm, "_", sp_nm, "_ensemble.grd")
        replace_grd_fx(file_name1)
        
        file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                            proj_nm, "_", sp_nm, "_ensemble_", eval_stat, "bin.grd")
        replace_grd_fx(file_name1)
        
        ####delete unused grids
        file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                            proj_nm, "_", sp_nm, "_", eval_stat, "bin.grd")
        unlink(file_name1, recursive=T, force=T)
        unlink(sub(pattern = ".grd$", replacement = ".gri", file_name1), recursive=T, force=T)
        
        file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                            proj_nm, "_", sp_nm, ".grd")
        unlink(file_name1, recursive=T, force=T)
        unlink(sub(pattern = ".grd$", replacement = ".gri", file_name1), recursive=T, force=T)
        
      }
    }
    
    #############################
    #delete temp raster files
    sp_nm = as.character(sp_nm) 
    sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
    temp_sp_files_to_delete<-paste0(project_path, sp_dir, "delete_temp_sp_files/", "*")
    unlink(temp_sp_files_to_delete, recursive=T, force=T) #delete previous frames
    #Loc <- "mydir"
    #system(paste0("rm -r ", temp_sp_files_to_delete))
    
    
    # return output to console
    sink(NULL) 
  }
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
##### END MODEL PROJECTIONS #####
#################################



temp_loc_to_delete=paste0("E:/Invasive_SDMs/global_model/temp/", "*")
unlink(temp_loc_to_delete, recursive=T, force=T) #delete previous frames
