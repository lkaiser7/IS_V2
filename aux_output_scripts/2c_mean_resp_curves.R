### mean ensemble modeling curves ###

# clear the environment, temp files, and all other variables
#rm(list = ls())

##################
##### SET UP #####
##################

### FROM MASTER SCRIPT ###

# set root path to source files
#rootDir<-"D:/Phase1_SDMs/"
# set working directory to main analysis folder
setwd(rootDir)

# list model scales
all_mod_scale = c("global_notHI", "local_HI", "nested_HI")#c("local", "global", "nested")
all_mod_scale=paste0(all_mod_scale, "_models")
# load necessary packages
library(biomod2)
library(stringr)
library(raster)
library(dplyr)
library(ggplot2)

# load code to override and set manual facet scales
source(paste0(codeDir, "aux_output_scripts/facet_scale_fix.R"))
# alternatively use cowplot or grid.arrange on multiple plots

# set all_sp_nm = 'Clidemia_hirta' for testing and debugging
# list all species names to be analyzed

##############################
##### LOAD AND LOOP DATA #####
##############################

# loop through each species
s=1
for(s in 1:length(all_sp_nm)){ # set s = 1 for debugging
  # select single species name
  sp.name<-all_sp_nm[s]
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp.name, "_", "."), "/")
  
  # print species name for counter
  print(sp.name)
  
  # loop through each model type 
  m=1
  for(m in 1:length(all_mod_scale)){ # set m = 1 for debugging 
    # set model scale
    mod_scale<-all_mod_scale[m]
    
    # set path of ongoing project run for all outputs
    project_path<-paste0(rootDir, mod_scale, "/")
    
    # location to save any extra or more specific outputs
    outDir<-paste0(project_path, "outputs/")
    
    # change working directory to project path to save model outputs
    setwd(project_path)
    
    # load variable importance file
    var_imp<-read.csv("outputs/all_VariImp_model_mean.csv")
    # select species model means
    sp_var_imp<-var_imp[which(var_imp$SPECIES == sp.name),]
    # add column for model type
    sp_var_imp$model_scale<-mod_scale
    
    # load model fitting data
    workspace_name = paste0(sp_dir, sp.name,"_modelfitting.RData")
    load(workspace_name)
    # myBiomodModelOut
    
    # load the models to extract the predicted response curves (220 runs)
    bm.mod.names<-BIOMOD_LoadModels(myBiomodModelOut)
    
    # get 2D response data (takes ~3 min per species for all 220 runs)
    rp.dat<-response.plot2(
      models = myBiomodModelOut@models.computed,
      Data = get_formal_data(myBiomodModelOut, 'expl.var'),
      show.variables = get_formal_data(myBiomodModelOut, 'expl.var.names'),
      data_species = get_formal_data(myBiomodModelOut, 'resp.var'),
      do.bivariate = FALSE,
      fixed.var.metric = 'median',
      col = c("blue", "red"),
      legend = FALSE,
      plot = FALSE)
    
    # View(rp.dat)
    
    # # select full models only
    # rp.full<-rp.dat[grep(x=rp.dat$pred.name, pattern = paste0("_Full")),]
    # rp.dat<-rp.full
    
    ################################
    ##### FORMAT DATA AND PLOT #####
    ################################
    
    # add model type column
    rp.dat$model_type=NA
    # list chosen models
    chosen_models=models_to_run #c('GBM', 'MAXENT.Phillips')
    
    # loop through and table by model
    for(chosen_model in chosen_models){
      rp.dat[grep(x=rp.dat$pred.name, pattern = paste0("_", chosen_model, "$")),
             "model_type"]=chosen_model
    }
    dput(names(rp.dat))
    rp.dat2=rp.dat[, c("expl.name", "expl.val", "pred.val", "model_type")]
    
    # summarize values by the new model_type expl.name and expl.val columns
    rp.dat2=aggregate(pred.val ~ model_type+expl.name+expl.val, data = rp.dat2, mean)
    # View(rp.dat2)
    
    # create folder in output directory for mean response curve plots
    rt_fold<-paste0(outDir, "mean_response_tables/")
    dir.create(rt_fold, showWarnings = FALSE)
    # save table for the global, local and nested run
    write.csv(rp.dat2, paste0(rt_fold, sp.name, "_pred_val_by_model.csv"))
    
    # save tables per model scale run
    if(mod_scale == "local_HI"){ #"global_notHI", "local_HI", "nested_HI"
      local_tab<-rp.dat2
    }else if(mod_scale == "global_notHI"){
      global_tab<-rp.dat2
    }else{
      nested_tab<-rp.dat2
    }
    
    # then merge to make a the panel figure
    
    # add model scale column
    rp.dat2$model_scale=mod_scale
    # merge data for plotting
    if(m == 1){
      all.rp.dat2<-rp.dat2
      all_var_imp<-sp_var_imp
    }else{
      all.rp.dat2<-rbind(all.rp.dat2, rp.dat2)
      all_var_imp<-rbind(all_var_imp, sp_var_imp)
    }
    
  } # END m loop
  
  # reset working directory to project
  setwd(rootDir)
  
  # create folder in output directory for mean response curve plots
  rc_fold<-paste0(rootDir, "combined_results/mean_response_curves/")
  dir.create(rc_fold, showWarnings = FALSE, recursive = T) 
  
  # color blind palettes
  # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible 
  cb_palette<-c("#009E73", "#E69F00", "#CC79a7") # "#56B4E9"
  
  # ORIGINAL CODE: this is not totally fixed, but nearly there
  # TO FIX, COLOR VARIABLE WILL BE MODEL_SCALE VARIABLE DENOTING GLOBAL, LOCAL, OR NESTED
  # ggplot(local_tab, 
  #        aes(x = expl.val, y = pred.val, colour = model_type, group = expl.name)) + 
  #   geom_line(size = 1) + facet_wrap(~ expl.name, scales = 'free_x') +
  #   labs(x = '', y = 'probability of occurrence', colour = 'model type') + 
  #   scale_color_brewer(type = 'qual', palette = 4) + ggtitle(sp.name) + theme_minimal() + 
  #   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  
  # rename variable importance column names
  names(all_var_imp)[3]<-"expl.name"
  
  ###################################
  # GBM Means
  gbm_data<-all.rp.dat2[which(all.rp.dat2$model_type == "GBM"),]
  # gbm_mean<-aggregate(pred.val ~ expl.name+model_scale, data = gbm_data, mean)
  # names(gbm_mean)[3]<-"pred.mean"
  # # View(gbm_data)
  
  # merge data for line size reference
  gbm_plot_data<-merge(gbm_data, all_var_imp, by =  c("expl.name", "model_scale"), all = T)
  gbm_plot_data$line_size<-gbm_plot_data$GBM_MEAN*3
  # dim(gbm_data); dim (gbm_plot_data); head(gbm_plot_data)
  
  ###################################
  # MaxEnt Means
  max_data<-all.rp.dat2[which(all.rp.dat2$model_type == "MAXENT.Phillips"),]
  # max_mean<-aggregate(pred.val ~ expl.name+model_scale, data = max_data, mean)
  # names(max_mean)[3]<-"pred.mean"
  # # max_mean
  
  # merge data for line size reference
  max_plot_data<-merge(max_data, all_var_imp, by =  c("expl.name", "model_scale"), all = T)
  max_plot_data$line_size<-max_plot_data$MAXENT_MEAN*3
  # dim(max_data); dim (max_plot_data); head(max_plot_data)
  
  
  ############
  #calculate squared deviation 
  #View(gbm_plot_data)
  #dput(names(gbm_plot_data))
  ooo=1
  for (ooo in c(1,2)){
    #View(tmp_df)
    mean_var_imp_df=gbm_plot_data[,c("expl.name", "model_scale", "MAXENT_MEAN", "GBM_MEAN")]
    library(reshape2)
    if (ooo==1){
      model_plot_data_reformat=gbm_plot_data[, c("SPECIES", "expl.name", "GBM_MEAN", "expl.val", "pred.val", 
                                                 "model_scale")]
      mean_var_imp_df=dcast(mean_var_imp_df, expl.name ~ model_scale, mean, value.var="GBM_MEAN")
      mod_type="GBM"
    }
    if (ooo ==2){
      model_plot_data_reformat=max_plot_data[, c("SPECIES", "expl.name", "GBM_MEAN", "expl.val", "pred.val", 
                                                 "model_scale")]
      mean_var_imp_df=dcast(mean_var_imp_df, expl.name ~ model_scale, mean, value.var="MAXENT_MEAN")
      
      mod_type="maxent"
    }
    names(mean_var_imp_df)=c("expl.name", "GvarImp", "LvarImp", "NvarImp")
    mean_var_imp_df$GL_meanImp=(mean_var_imp_df$GvarImp+mean_var_imp_df$LvarImp)/2
    mean_var_imp_df$GN_meanImp=(mean_var_imp_df$GvarImp+mean_var_imp_df$NvarImp)/2
    mean_var_imp_df$LN_meanImp=(mean_var_imp_df$LvarImp+mean_var_imp_df$NvarImp)/2
    #first step is to standardize the expl var values as the response curves do not yield equal values!
    bio = var_names[1]
    for (bio in var_names){
      model_plot_data_reformat_single=model_plot_data_reformat[model_plot_data_reformat$expl.name==bio,]
      #View(model_plot_data_reformat_single)
      min_expl_val=max(as.matrix(dcast(model_plot_data_reformat_single, formula = . ~ model_scale, min, value.var="expl.val")[,-1]))
      max_expl_val=min(as.matrix(dcast(model_plot_data_reformat_single, formula = . ~ model_scale, max, value.var="expl.val")[,-1]))
      #create a sequence only for values common across all 3 scales
      val_seq=seq(min_expl_val, max_expl_val, length.out=100)
      
      model_scale = all_mod_scale[1]
      for (model_scale in all_mod_scale){
        model_plot_data_reformat_single_bio_and_scale=model_plot_data_reformat_single[model_plot_data_reformat_single$model_scale==model_scale,]
        extrapolated_response=approx(model_plot_data_reformat_single_bio_and_scale$expl.val, model_plot_data_reformat_single_bio_and_scale$pred.val, xout=val_seq)
        if (model_scale == all_mod_scale[1]){
          extrapolated_response_DF=data.frame(expl=val_seq, model_scale=extrapolated_response$y)
        }else{
          extrapolated_response_DF=cbind(extrapolated_response_DF, model_scale=extrapolated_response$y)
        }
      }
      names(extrapolated_response_DF)=c("expl_vals", all_mod_scale)
      extrapolated_response_DF$GL_diff=abs(extrapolated_response_DF$global_notHI_models-extrapolated_response_DF$local_HI_models)
      extrapolated_response_DF$GN_diff=abs(extrapolated_response_DF$global_notHI_models-extrapolated_response_DF$nested_HI_models)
      extrapolated_response_DF$LN_diff=abs(extrapolated_response_DF$local_HI_models-extrapolated_response_DF$nested_HI_models)
      #View(extrapolated_response_DF)
      extrapolated_response_DF_tmp=extrapolated_response_DF[,c("expl_vals", "global_notHI_models", "local_HI_models", "nested_HI_models")]
      extrapolated_response_DF_tmp=melt(extrapolated_response_DF_tmp,id.vars = "expl_vals")
      ggplot(extrapolated_response_DF_tmp, aes(x = expl_vals, y = value, colour = variable))+geom_line()
      
      extrapolated_response_DF_tmp=extrapolated_response_DF[,c("expl_vals", "GL_diff", "GN_diff", "LN_diff")]
      extrapolated_response_DF_tmp=melt(extrapolated_response_DF_tmp,id.vars = "expl_vals")
      ggplot(extrapolated_response_DF_tmp, aes(x = expl_vals, y = value, colour = variable))+geom_line()
  
      # SS_G_L_diff=sum(extrapolated_response_DF$G_L_diff^2)
      # SS_G_N_diff=sum(extrapolated_response_DF$G_N_diff^2)
      # SS_L_N_diff=sum(extrapolated_response_DF$L_N_diff^2)
      
      SS_GL_diff=sum(extrapolated_response_DF$GL_diff)
      SS_GN_diff=sum(extrapolated_response_DF$GN_diff)
      SS_LN_diff=sum(extrapolated_response_DF$LN_diff)
      
      bio_mean_var_imp_df=mean_var_imp_df[mean_var_imp_df$expl.name==bio,-1]
      results_row=data.frame(model=mod_type, species=sp.name, bio, SS_GL_diff, SS_GN_diff, SS_LN_diff) 
      results_row=cbind(results_row, bio_mean_var_imp_df)
      
      if (ooo==1 & s==1 & bio==var_names[1]){
        SS_results_DF=results_row
      }else{
        SS_results_DF=rbind(SS_results_DF, results_row)
      }
    }
  }
  #View(SS_results_DF)
  ggplot(gbm_plot_data,
         aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
    geom_line(lwd = gbm_plot_data$line_size) + #geom_line(lwd = gbm_plot_data$pred.mean) + 
    facet_wrap(~ expl.name, scales = 'free_x') +
    # coord_cartesian(xlim = c(0, 150)) +  # sets manual x scale limits
    facet_wrap_custom(~ expl.name, scales = "free", scale_overrides = list(
      scale_override(1, scale_x_continuous(limits = c(0, 50))),
      scale_override(3, scale_x_continuous(limits = c(0, 150))),
      scale_override(4, scale_x_continuous(limits = c(0, 40)))
    )) +
    labs(x = '', y = 'probability of occurrence', colour = 'model scale') + 
    # scale_color_brewer(type = 'qual', palette = 4) + 
    scale_color_manual(values = cb_palette) + 
    ggtitle(paste(sub("_", " ", sp.name), "GBM")) + 
    theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(rc_fold, sp.name, "_GBM_mean_resp_curve.tiff"))
  
  # MAXENT Plot
  ggplot(max_plot_data,
         aes(x = expl.val, y = pred.val, colour = model_scale), group = expl.name) + 
    geom_line(lwd = max_plot_data$line_size) + #geom_line(lwd = max_plot_data$pred.mean) + 
    facet_wrap(~ expl.name, scales = 'free_x') +
    # coord_cartesian(xlim = c(0, 150)) +  # sets manual x scale limits
    facet_wrap_custom(~ expl.name, scales = "free", scale_overrides = list(
      scale_override(1, scale_x_continuous(limits = c(0, 50))),
      scale_override(3, scale_x_continuous(limits = c(0, 150))),
      scale_override(4, scale_x_continuous(limits = c(0, 40)))
    )) +
    labs(x = '', y = 'probability of occurrence', colour = 'model scale') + 
    # scale_color_brewer(type = 'qual', palette = 4) + 
    scale_color_manual(values = cb_palette) + 
    ggtitle(paste(sub("_", " ", sp.name), "MaxEnt")) + 
    theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(rc_fold, sp.name, "_MAXENT_mean_resp_curve.tiff"))
  
} # END s loop
write.csv(SS_results_DF, paste0(rc_fold,"_response_deviations.csv"), row.names = F)



##### END MEAN RESPONSE CURVES #####
