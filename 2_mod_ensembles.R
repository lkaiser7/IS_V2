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

# reset species counter to 1
sp_nm = all_sp_nm[1] 

# initialize snowfall parallel computing function
sp_parallel_run = function(sp_nm){  
  # load necessary packages
  library(biomod2)
  library(stringr)
  require(snowfall)
  library(raster)
  
  # convert species name to character object (in case the species are numbered)
  sp_nm = as.character(sp_nm) 
  sp_nm0=sp_nm
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
  #gc()
  # print posting of temporary file location
  cat('\n temporary files to be deleted saved here:', temp_sp_files_to_delete, '\n')
  
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
    if (packageVersion("biomod2")=="4.0"){
      myBiomodEM<-BIOMOD_EnsembleModeling( 
        bm.mod = myBiomodModelOut,  #BIOMOD.models.out from model fitting
        models.chosen = remaining_models,  #vector of model runs to use
        em.by = 'all',  #how models will be combined
        metric.eval = eval_stats, #evaluation metrics to build ensemble
        metric.select.thresh = eval.metric.threshold,  #threshold to exclude models
        prob.mean = TRUE,  #estimate mean probabilities 
        prob.cv = TRUE,  #estimate coefficient of variation
        prob.ci = TRUE,  #estimate confidence interval of prob.mean
        prob.ci.alpha = 0.05,  #signficance level for estimating confidence interval
        prob.median = TRUE,  #estimate median
        committee.averaging = TRUE,  #estimate committee averaging
        prob.mean.weight = TRUE,  #estimate weighted sums
        prob.mean.weight.decay = 'proportional' ) #define relative importance of weights
    }else{
      myBiomodEM<-BIOMOD_EnsembleModeling( 
        modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
        chosen.models = remaining_models, #"all", #remaining_models,  #vector of model runs to use
        em.by = 'all', # 'all' 'PA_dataset+repet', 'PA_dataset',   #how models will be combined into ensemble. If all, will create a single ensemble model https://rpubs.com/dgeorges/38564 #if all, evaluation is using all candidate PAs
        #em.by = 'PA_dataset', # 'all' 'PA_dataset+repet', 'PA_dataset', 'PA_dataset+algo'#'all',  #how models will be combined #CHECK THIS: https://rpubs.com/dgeorges/38564 #if all, evaluation is using all candidate PAs
        eval.metric = eval_stats, #evaluation metrics to build ensemble
        eval.metric.quality.threshold = NULL, #eval.metric.threshold,  #threshold to exclude models
        prob.mean = T,  #estimate mean probabilities 
        prob.cv = F,  #estimate coefficient of variation
        prob.ci = F,  #estimate confidence interval of prob.mean
        prob.ci.alpha = 0.05,  #signficance level for estimating confidence interval
        prob.median = F,  #estimate median
        committee.averaging = T,  #estimate committee averaging
        prob.mean.weight = TRUE,  #estimate weighted sums
        prob.mean.weight.decay = 'proportional' ) #define relative importance of weights
    }
    
    # reset working directory to root 
    # setwd(rootDir)
    
    # sign-posting of completed ensemble model creation
    cat('\n', sp_nm, 'ensemble complete.') 
    # record time and date stamp
    cat(format(Sys.time(), "%a %b %d %Y %X"))
    
    ########################
    ### ENSEMBLE OUTPUTS ###
    ########################
    # print ensemble modeling summary
    cat('\n', sp_nm, ': start ensemble outputs') 
    
    myBiomodEM
    
    ##################
    ##################
    cat('\n', sp_nm, ': create varImp for entire ensemble') 
    cat('\n', sp_nm, ': get ensemble model name') 
    
    # need to complete this below
    # return output model evaluation metrics results
    myBiomodEMEval<-get_evaluations(myBiomodEM)
    myBiomodEMEval=myBiomodEMEval[[paste0(sp_nm, "_EM",spp_ensemble_type,"By",eval_stats,"_mergedAlgo_mergedRun_mergedData")]]
    myBiomodEMEval=myBiomodEMEval[row.names(myBiomodEMEval)==eval_stats,]
    # assign file path name for results
    FileName<-paste0(project_path, sp_dir, sp_nm, "_",eval_stats,"_EM.csv")
    # create .csv file and save TSS outputs
    write.table(myBiomodEMEval, file = FileName, sep = ",", col.names = NA)
    
    cat('\n', sp_nm, ': calculate variable importance') 
    
    load(paste0(project_path, sp_dir, sp_nm0, "_BiomodData.RData")) #myBiomodData
    
    EM_model=myBiomodEM@em.models[[paste0(sp_nm, "_EM",spp_ensemble_type,"By",eval_stats,"_mergedAlgo_mergedRun_mergedData")]]
    # change working directory to project path to save model outputs
    #setwd(project_path)
    EM_var_imp=variables_importance(model = EM_model, data = myBiomodData@data.env.var, nb_rand =2)
    EM_var_imp=EM_var_imp$mat
    EM_var_imp=data.frame(pred=row.names(EM_var_imp), varImp=apply(EM_var_imp, 1, mean))
    #names(EM_var_imp)=c("pred", "varImp")
    #setwd(rootDir)
    
    FileName00=paste0(project_path, sp_dir, sp_nm, "_",eval_stats,"_EMVarImp.csv")
    write.table(EM_var_imp, file = FileName00, sep = ",", col.names = NA)
    
    #############################
    cat('\n', sp_nm, ': create ensemble response curve') 
    #setwd(project_path)
    curve_models=BIOMOD_LoadModels(myBiomodEM)
    curve_models=curve_models[grep(pattern = paste0(spp_ensemble_type,"By", eval_stats), x = curve_models)]
    #curve_models
    
    tiff_name=paste0(project_path, sp_dir, sp_nm, "_",eval_stats,"_EM_response_curveV2.tif")
    

    create_response_plots_fx=function(regression_mat=myBiomodData@data.env.var, fit=get(curve_models[1]), tiff_name){
      #response curve data using transformed matrix
      resp_curve_data0=regression_mat
      df_classes=sapply(resp_curve_data0, class)
      numeric_cols=names(resp_curve_data0)[df_classes %in% c("numeric", "integer")]
      
      numeric_col=numeric_cols[1]
      unique_factor_val=1
      for (numeric_col in numeric_cols){
        cat("doing ", numeric_col, "\n")
        
        #create median matrix
        jnk=as.data.frame(t(as.matrix(apply(resp_curve_data0[,c(numeric_cols)], 2, FUN = median))))
        resp_curve_data_frst_row=resp_curve_data0[1,]
        col_match=match(names(jnk), names(resp_curve_data_frst_row))
        resp_curve_data_frst_row[,col_match]=jnk
        resp_curve_data=resp_curve_data_frst_row[rep(1, times=100),]
        
        #now add the sequence of values from the target response variable
        #use the untransformed sequence!
        response_var_Q99_notTrans=quantile(resp_curve_data0[,numeric_col], 0.99)
        response_var_Q01_notTrans=quantile(resp_curve_data0[,numeric_col], 0.01)
        resp_var_seq=seq(from=min(resp_curve_data0[,numeric_col],na.rm=T), to=max(resp_curve_data0[,numeric_col],na.rm=T), length.out=100)
        resp_curve_data[,numeric_col]=resp_var_seq
        
        #now the df to use for prediction is ready
        resp_curve_vals=predict(fit, newdata = resp_curve_data)
        resp_curve_results=cbind(resp_curve_data, response=resp_curve_vals)
        resp_curve_results=resp_curve_results[,c(numeric_col, "response")]
        
        resp_curve_results = aggregate(resp_curve_results,
                                       by = list(resp_curve_results[, numeric_col]),
                                       FUN = mean)
        resp_curve_results=resp_curve_results[, c(numeric_col, "response")]
        resp_curve_results$pred_name=numeric_col
        names(resp_curve_results)[1]="pred"
        
        if (numeric_col == numeric_cols[1]){
          all_resp_curve_results=resp_curve_results
        }else{
          all_resp_curve_results=rbind(all_resp_curve_results, resp_curve_results)
        }
      }
      
      tiff(tiff_name,
           width = 8, height = 8, units = "in",
           pointsize = 12, compress="lzw", bg = "white", res = 300)
      par(mfrow=c(2,2))
      for (numeric_col in numeric_cols){
        tmp_all_resp_curve_results=all_resp_curve_results[all_resp_curve_results$pred_name==numeric_col,]
        plot(tmp_all_resp_curve_results$pred, tmp_all_resp_curve_results$response, xlab=numeric_col, ylab="", ylim=c(0,1), type="l", col="darkblue", lwd=2)  
      }
      graphics.off()
      return(all_resp_curve_results)
    }
    # graphics.off()
    response_plot_data=create_response_plots_fx(regression_mat=myBiomodData@data.env.var, fit=get(curve_models[1]), tiff_name)
    
    cat('\n', sp_nm, ': write response curve data') 
    #response_plot_data=response_plot_data[,-c(1, 4)]
    #View(response_plot_data)
    FileName00=paste0(project_path, sp_dir, sp_nm, "_",eval_stats,"_EM_response_curve.csv")
    write.csv(response_plot_data, FileName00, row.names = F)
    #setwd(rootDir)
    
    ##################
    ##################
    #reset working directory to root 
    setwd(rootDir)
    
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
    # reset sink to console output
    sink(NULL)
  }else{
    # sign-posting if file for ensemble modeling has already been created 
    cat('\n', sp_nm, 'ensemble modleling already done...')
    # indicates that this species has already been run 
    # reset sink to console output
    sink(NULL)
    unlink(txt_nm) #delete previous frames
  }
  
  # delete select temporary files per species once processing is finished
  unlink(paste0(temp_sp_files_to_delete, "*"), recursive=T, force=T) #delete previous frames
  file.remove(list.files(tempdir(), full.names = T, pattern = "^file")) #https://stackoverflow.com/questions/45894133/deleting-tmp-files
  gc()
  
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
system.time((sfLapply(all_sp_nm, fun = sp_parallel_run)))
# remove all global environmental variables from cluster cores
sfRemoveAll()
# end parallel computing on other cpu cores
sfStop()


#################################
##### END ENSEMBLE MODELING #####
#################################

