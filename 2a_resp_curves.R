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

# loop through all species
sp_nm = all_sp_nm[1]
for (sp_nm in all_sp_nm){
  # convert species names to characters
  sp_nm = as.character(sp_nm)  
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm, "_", "."), "/")
  
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
    tiff_name=paste0(rc_fold, sp_nm, "_response_curve_", model, "all_vars.tif")
    # create blank image file
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
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
      dev.off()
      ###These next sections need to be changed as "myRespPlot2D" is now a list###
      
      # # set for debugging
      # bioclim_cnt = 1
      # preds=unique(myRespPlot2D$expl.name)
      # # loop through each bioclimatic variable
      # pred = preds[1]
      # for (pred in preds) {
      #   pred_myRespPlot2D=myRespPlot2D[myRespPlot2D$expl.name==pred,]
      #   # maximum y limit
      #   ymax_lim = 0
      #   # minimum y limit
      #   ymin_lim = 10000
      #   # loop through all repeated model runs
      #   for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
      #     # maximums of y limit per variable
      #     ymax_lim = max(ymax_lim, max(myRespPlot2D[[bioclim_cnt]][rep]))
      #     # minimums of y limit per variable
      #     ymin_lim = min(ymin_lim, min(myRespPlot2D[[bioclim_cnt]][rep]))
      #   }
      #   # maximum x limit
      #   xmax_lim = max(myRespPlot2D[[bioclim_cnt]][1])
      #   # minimum x limit
      #   xmin_lim = min(myRespPlot2D[[bioclim_cnt]][1])
      #   # list variable name
      #   var_name = names(myRespPlot2D)[bioclim_cnt]
      #   
      #   # name output image file
      #   tiff_name=paste0(rc_fold, sp_nm, "_response_curve_", model, "_", var_name, ".tif")
      #   # create blank image file for each bioclimatic variable
      #   tiff(tiff_name, res = 300, units = "in", pointsize = 12,
      #        width = 10, height = 10, compression = "lzw")
      #   
      #   # get mean of 2D response plot per variable
      #   var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
      #   # loop through each model run per variable
      #   for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
      #     # get mean of all other model runs and variables
      #     pred = rowMeans(myRespPlot2D[[bioclim_cnt]][rep])
      #     # for following plots after first model run
      #     if (rep == 2){
      #       # plot variables per prediction
      #       plot(var, pred, type = "l", col = "grey", 
      #            xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), 
      #            xlab = var_name, ylab = "Response")
      #     } else {
      #       # add lines
      #       lines(var, pred, type = "l",  col = "grey", 
      #             xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim))
      #     }
      #   }
      #   
      #   # get mean of 2D response plot per variable
      #   var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
      #   # get mean of all other 2D response plots and variables
      #   pred = rowMeans(myRespPlot2D[[bioclim_cnt]][2:length(myRespPlot2D[[bioclim_cnt]])])
      #   # add average response lines
      #   lines(var, pred, type = "l", lwd = 3, 
      #         xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim))    
      #   # save image file
      #   dev.off()
      # }
      # 
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
      dev.off()
      # # load and format data as done previously done above for single variable plots
      # ymax_lim = 0
      # ymin_lim = 10000
      # for (bioclim_cnt in 1:length(myRespPlot2D)) {
      #   for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
      #     ymax_lim = max(ymax_lim, max(myRespPlot2D[[bioclim_cnt]][rep]))
      #     ymin_lim = min(ymin_lim, min(myRespPlot2D[[bioclim_cnt]][rep]))
      #   }
      # }
      # 
      # # set graphical parameters for comination plot
      # par(mfrow = c(2,2), oma = c(0,0,3,0))
      # # load and format data as done previously done above to create combo plot
      # bioclim_cnt = 1
      # for (bioclim_cnt in 1:length(myRespPlot2D)) {
      #   xmax_lim = max(myRespPlot2D[[bioclim_cnt]][1])
      #   xmin_lim = min(myRespPlot2D[[bioclim_cnt]][1])
      #   var_name = names(myRespPlot2D)[bioclim_cnt]
      #         
      #   for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
      #     var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
      #     pred = rowMeans(myRespPlot2D[[bioclim_cnt]][rep])
      #     if (rep == 2){
      #       plot(var, pred, type = "l", , col = "grey",
      #            xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), 
      #            xlab = var_name, ylab = "Response")
      #     } else {
      #       lines(var, pred, type = "l", col = "grey", 
      #             xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim))
      #     }
      #   }
      #   
      #   var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
      #   pred = rowMeans(myRespPlot2D[[bioclim_cnt]][2:length(myRespPlot2D[[bioclim_cnt]])])
      #   lines(var, pred, type="l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), lwd = 3)          
      # } 
      # # add title for combo plot
      # temp_title = paste0(model, " modeled response for ", sp_nm)
      # # place title on plot
      # mtext(temp_title, adj = 0.5, side = 3, outer = TRUE) 
      # save image file
      
      # sign-posting of completed response curve plotting
      cat('\n response curves for ',sp_nm, model, 'finished.')  
    }else{
      # sign-posting if file for response curves has already been created 
      cat('\n response curves for ',sp_nm, model, 'already done.')  
    }      
  }
}

# reset working directory to root 
setwd(rootDir)

####################################
##### END RESPONSE CURVE PLOTS #####
####################################