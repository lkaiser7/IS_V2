# load necessary packages
require(snowfall)

# reset sp_nm counter to first species
sp_nm = all_sp_nm[1] 

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
  
  # convert species name to character object (in case the species are numbered)
  sp_nm = as.character(sp_nm) 
  sp.nm=sub("_", ".", sp_nm)
  
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  
  # create temporary directory per species to be deleted
  temp_sp_files_to_delete<-paste0(project_path, sp_dir, "delete_temp_sp_files/")
  dir.create(temp_sp_files_to_delete, showWarnings = FALSE)
  # set temporary directory to created temp file
  rasterOptions(tmpdir = temp_sp_files_to_delete)
  unlink(paste0(temp_sp_files_to_delete, "*"), recursive=T, force=T) #delete previous temp files if any
  file.remove(list.files(tempdir(), full.names = T, pattern = "^file")) #https://stackoverflow.com/questions/45894133/deleting-tmp-files
  #gc()
  # print posting of temporary file location
  cat('\n temporary files to be deleted saved here:', temp_sp_files_to_delete, '\n')
  
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
  workspace_name_out = paste0(project_path, sp_dir, sp_nm, "_", proj_nm, "_projections.RData") 
  
  #####################
  ##### LOAD DATA #####
  #####################
  
  # sign-posting of loading bioclimatic variables data
  cat('\n loading predictors...')    
  # sign-posting indicating which raster is used for projections
  cat('\n using these env files for projection raster:', env_var_files, 
      '\n from dir:', clim_data)
  
  eval_projects=c("regional_HI", "global_notHI", "nested_HI")
  eval_project=eval_projects[1]
  for (eval_project in eval_projects){
    n_accuraccy_iters=10
    if (eval_project=="nested_HI"){
      reg_project_path=paste0(rootDir,"regional_HI","_models/")
      projection_name = paste0(eval_project, "_accuracy")
      predictors=read.csv(paste0(reg_project_path, sp_dir, sp_nm, "_bioclim_points.csv"))
      #get equal amounts of presence/ absence data
      tmp_pres=predictors[predictors$PA %in% 1,]
      n_pres=nrow(tmp_pres)
      tmp_pas=predictors[!(predictors$PA %in% 1),]
      tmp_pas=tmp_pas[sample(nrow(tmp_pas), n_pres*n_accuraccy_iters),]
      reg_predictors=rbind(tmp_pres, tmp_pas)
      
      reg_project_path=paste0(rootDir,"global_notHI","_models/")
      projection_name = paste0(eval_project, "_accuracy")
      predictors=read.csv(paste0(reg_project_path, sp_dir, sp_nm, "_bioclim_points.csv"))
      #get equal amounts of presence/ absence data
      tmp_pres=predictors[predictors$PA %in% 1,]
      n_pres=nrow(tmp_pres)
      tmp_pas=predictors[!(predictors$PA %in% 1),]
      tmp_pas=tmp_pas[sample(nrow(tmp_pas), n_pres*n_accuraccy_iters),]
      glo_predictors=rbind(tmp_pres, tmp_pas)
      
      predictors=rbind(reg_predictors, glo_predictors)
    }else{
      reg_project_path=paste0(rootDir,eval_project,"_models/")
      projection_name = paste0(eval_project, "_accuracy")
      predictors=read.csv(paste0(reg_project_path, sp_dir, sp_nm, "_bioclim_points.csv"))
      #get equal amounts of presence/ absence data
      tmp_pres=predictors[predictors$PA %in% 1,]
      n_pres=nrow(tmp_pres)
      tmp_pas=predictors[!(predictors$PA %in% 1),]
      tmp_pas=tmp_pas[sample(nrow(tmp_pas), n_pres*n_accuraccy_iters),]
      predictors=rbind(tmp_pres, tmp_pas)
    }
    
    actual=predictors$PA
    predictors=predictors[,var_names]
    
    # sign-posting noting completion of loading predictors
    cat('\n done loading predictors.')  
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
    
    # use previously created baseline projections
    myBiomodProjection<-myBiomodProj_baseline
    
    # sign-posting to run forecasting of ensemble models
    cat('\n run ensemble forecasting...')
    ################################
    ##### ensemble_forecasting #####
    ################################
    
    # change working directory to project path to save model outputs
    setwd(project_path)
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      projection.output = myBiomodProjection,  #BIOMOD.projection.out from projections
      total.consensus = TRUE,  #mean of all combined model projections
      selected.models= "all", #remaining_models,
      #new.env = predictors,
      EM.output = myBiomodEM, #BIOMOD.EnsembleModeling.out from ensemble modeling
      proj.name=proj_nm,
      binary.meth = eval_stats,  #evaluation method statistics 
      keep.in.memory = memory) #,
    
    EM_projections=as.data.frame(myBiomodEF@proj@val)
    EM_projections=EM_projections[,grep(pattern = spp_ensemble_type, names(EM_projections))] 
    EM_projections=data.frame(actual,EM_projections)
    EM_projections[is.na(EM_projections$actual),"actual"]=0
    #make equal number of pseudo abs
    # n_pres=sum(EM_projections$actual==1)
    # tmp_pres=EM_projections[EM_projections$actual==1,]
    # tmp_pas=EM_projections[EM_projections$actual==0,]
    # tmp_pas=tmp_pas[sample(nrow(tmp_pas), n_pres),]
    # EM_projections=rbind(tmp_pres, tmp_pas)
    #View(EM_projections)
    
    #################################
    #calculate own cutoff!
    eval_stat = eval_stats[1]
    for (eval_stat in eval_stats){
      scores_all=get_evaluations(myBiomodEM)
      cat("doing ", eval_stat, " ensemble bin rasters \n")
      index=grep(pattern = paste0("EMwmeanBy", eval_stat),names(scores_all))
      jnk=scores_all[[index]]
      cutoff=jnk[eval_stat, "Cutoff"]
      iter=1
      for (iter in c(1:n_accuraccy_iters)){
        EM_projections_partition=EM_projections[EM_projections$actual==1,]
        tmp=EM_projections[c((iter*n_pres+1):(iter*n_pres+n_pres)),]
        EM_projections_partition=rbind(EM_projections_partition, tmp)
        #View(EM_projections_partition)
        custom_cutoff=biomod2:::Find.Optim.Stat(Stat = "TSS", Fit=EM_projections_partition$EM_projections, Obs=EM_projections_partition$actual, )[2]
        #################################
        #binarize
        #get threshold values and apply to appropriate ensemble layer
        cut_off_to_use=custom_cutoff
        EM_projections_bin=EM_projections_partition
        EM_projections_bin$EM_projections=as.numeric(EM_projections_bin$EM_projections>=cut_off_to_use)
        conf_matrix<-table(EM_projections_bin$EM_projections,EM_projections_bin$actual)
        library(caret)
        sp_sens=sensitivity(conf_matrix)
        sp_spec=specificity(conf_matrix)
        conf_mat_stats=confusionMatrix(conf_matrix)
        sp_accuracy=conf_mat_stats$overall[1]
        #View(EM_projections_bin)
        sp_accuracy_df=data.frame(eval_project, eval_stat, sp_nm, sp_accuracy, sp_sens, sp_spec, cutoff, custom_cutoff, cut_off_to_use)
        if (iter==1){
          custom_cutoffs=custom_cutoff
          spp_accuracy_df=sp_accuracy_df
        }else{
          custom_cutoffs=c(custom_cutoffs, custom_cutoff)
          spp_accuracy_df=rbind(spp_accuracy_df, sp_accuracy_df)
        }
      }
      sp_accuracy_df_pt1= spp_accuracy_df[1,c("eval_project", "eval_stat", "sp_nm")]
      row.names(sp_accuracy_df_pt1)=NULL
      sp_accuracy_df_pt2= spp_accuracy_df[,c("sp_accuracy", "sp_sens","sp_spec", "cutoff", "custom_cutoff", "cut_off_to_use")]
      sp_accuracy_df_pt2=apply(sp_accuracy_df_pt2,2,mean, na.rm=T)
      sp_accuracy_df_pt2=t(as.data.frame(sp_accuracy_df_pt2))
      row.names(sp_accuracy_df_pt2)=NULL
      sp_accuracy_df=cbind(sp_accuracy_df_pt1, sp_accuracy_df_pt2)
      FileName00=paste0(project_path, sp_dir, eval_project, "_", sp_nm, "_",eval_stat,"_EM_accuracy_DF.csv")
      write.csv(sp_accuracy_df, FileName00, row.names = F)
      
      # reset working directory to root 
      
    }
    setwd(rootDir)
  }
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

model_scales=c("global_notHI", "regional_HI", "nested_HI")
# model_scales=c("regional_HI", "nested_HI")
model_scale=model_scales[1]
for (model_scale in model_scales){
  project_run<-paste0(model_scale, "_models")
  # set path of ongoing project run for all outputs
  project_path<-paste0(rootDir, project_run, "/")
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
  
}


#################################
##### END MODEL PROJECTIONS #####
#################################

model_scale=model_scales[1]
for (model_scale in model_scales){
  tmp_project_path=paste0(rootDir, model_scale, "_models/")
  cat("doing ", model_scale, "\n")
  sp_nm = all_sp_nm[1]
  for (sp_nm in all_sp_nm){
    cat("doing ", sp_nm, "\n")
    eval_projects=c("regional_HI", "global_notHI", "nested_HI")
    eval_project=eval_projects[2]
    for (eval_project in eval_projects){
      eval_stat = eval_stats[1]
      for (eval_stat in eval_stats){
        # replace species naming convention of "_" with "." 
        sp.name_period=str_replace_all(sp_nm, "_", ".")
        sp_dir = paste0(sp.name_period, "/")
        
        FileName00=paste0(tmp_project_path, sp_dir, eval_project, "_", sp_nm, "_",eval_stat,"_EM_accuracy_DF.csv")
        sp_evalDF=read.csv(FileName00)
        sp_evalDF=cbind(model_scale, sp_evalDF)
        if (sp_nm == all_sp_nm[1] & model_scale==model_scales[1] & 
            eval_stat==eval_stats[1] & eval_project==eval_projects[1]){
          all_sp_evalDF=sp_evalDF
        }else{
          all_sp_evalDF=rbind(all_sp_evalDF, sp_evalDF)
        }
      }  
    }
  }  
}
all_sp_evalDF=cbind(all_sp_evalDF, current_spp_name=replace_spp_names(all_sp_evalDF$sp_nm))
#View(all_sp_evalDF)
write.csv(all_sp_evalDF, paste0(rootDir, "combined_results/all_reg_and_global_EM_eval.csv"))
