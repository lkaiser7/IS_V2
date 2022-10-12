### mean ensemble modeling curves ###
##################
##### SET UP #####
##################

### FROM MASTER SCRIPT ###
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

######################
#first collect EM var imp
sp.name = all_sp_nm[1]
model_scales=c("global_notHI", "local_HI", "nested_HI")
model_scale=model_scales[1]
for (model_scale in model_scales){
  tmp_project_path=paste0(rootDir, model_scale, "_models/")
  cat("doing ", model_scale, "\n")
  for (sp_nm in all_sp_nm){
    cat("doing ", sp_nm, "\n")
    
    # replace species naming convention of "_" with "." 
    sp.name_period=str_replace_all(sp_nm, "_", ".")
    sp_dir = paste0(sp.name_period, "/")
    
    FileName00=paste0(tmp_project_path, sp_dir, sp.name_period, "_",eval_stats,"_EMVarImp.csv")
    sp_varImpDF=read.csv(FileName00)
    sp_varImpDF=sp_varImpDF[,-1]
    sp_varImpDF$scale=model_scale
    sp_varImpDF$species=sp_nm
    if (sp_nm == all_sp_nm[1] & model_scale==model_scales[1]){
      all_sp_varImpDF=sp_varImpDF
    }else{
      all_sp_varImpDF=rbind(all_sp_varImpDF, sp_varImpDF)
    }
  }  
  #View(all_sp_varImpDF)
}
write.csv(all_sp_varImpDF, paste0(rootDir, "combined_results/all_EM_var_imp.csv"))

# create folder in output directory for mean response curve plots
rc_fold<-paste0(rootDir, "combined_results/EM_mean_response_curves/")
dir.create(rc_fold, showWarnings = FALSE, recursive = T) 

emVI_fold<-paste0(rootDir, "combined_results/EM_mean_VariImp_plots/")
dir.create(emVI_fold, showWarnings = FALSE, recursive = T) 

# loop through each species
s=1
for(s in 1:length(all_sp_nm)){ # set s = 1 for debugging
  # select single species name
  sp.name<-all_sp_nm[s]
  # replace species naming convention of "_" with "." 
  sp.name_period=str_replace_all(sp.name, "_", ".")
  sp_dir = paste0(sp.name_period, "/")
  
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
    
    rp.dat2=read.csv(paste0(project_path, sp_dir, sp.name_period, "_TSS_EM_response_curve.csv"))
    rp.dat2=rp.dat2[,c("pred_name", "pred", "response")]
    names(rp.dat2)=c("expl.name", "expl.val", "pred.val")
    #View(rp.dat2)
    
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
      #all_var_imp<-sp_var_imp
    }else{
      all.rp.dat2<-rbind(all.rp.dat2, rp.dat2)
      #all_var_imp<-rbind(all_var_imp, sp_var_imp)
    }
    
  } # END m loop
  
  # color blind palettes
  # https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible 
  cb_palette<-c("#009E73", "#E69F00", "#CC79a7") # "#56B4E9"
  
  #all.rp.dat2$SPECIES=sp.name
  #process varImp
  sp_varImpDF=all_sp_varImpDF[,c("species", "pred", "varImp", "scale")]
  sp_varImpDF=sp_varImpDF[sp_varImpDF$species==sp.name,]
  names(sp_varImpDF)=c("SPECIES", "expl.name", "varImp", "model_scale")
  sp_varImpDF=sp_varImpDF[,-1]
  sp_varImpDF$model_scale=paste0(sp_varImpDF$model_scale, "_models")
  plot_data<-merge(all.rp.dat2, sp_varImpDF, by =  c("expl.name", "model_scale"), all = T)
  plot_data$line_size<-plot_data$varImp*3
  #View(plot_data)
  
  ############
  #calculate squared deviation 
  #View(gbm_plot_data)
  #View(tmp_df)
  mean_var_imp_df=sp_varImpDF[,c("expl.name", "model_scale", "varImp")]
  
  library(reshape2)
  model_plot_data_reformat=plot_data[, c("expl.name", "varImp", "expl.val", "pred.val", 
                                         "model_scale")]
  mean_var_imp_df=dcast(mean_var_imp_df, expl.name ~ model_scale, mean, value.var="varImp")
  mod_type="EM"
  
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
    
    if (s==1 & bio==var_names[1]){
      SS_results_DF=results_row
    }else{
      SS_results_DF=rbind(SS_results_DF, results_row)
    }
  }

  #View(SS_results_DF)
  ggplot(plot_data,
         aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
    geom_line(lwd = plot_data$line_size) + #geom_line(lwd = plot_data$pred.mean) + 
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
    ggtitle(paste(sub("_", " ", sp.name))) + 
    theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve.tiff"))
  
  ##############################
  #no nested:
  #View(plot_data)
  plot_data_no_nested=plot_data[plot_data$model_scale!="nested_HI_models",]
  
  ggplot(plot_data_no_nested,
         aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
    geom_line(lwd = plot_data_no_nested$line_size) + #geom_line(lwd = plot_data$pred.mean) + 
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
    ggtitle(paste(sub("_", " ", sp.name))) + 
    theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve_GL.tiff"))
  
  
  ########
  #revised var imp plots
  # color blind palette
  cb_palette<-c("#009E73", "#E69F00", "#CC79a7", "#56B4E9")
  sp_varImpDF_for_graph=sp_varImpDF
  
  to_replace=c("global_notHI_models", "local_HI_models", "nested_HI_models")
  replacement=c("Global", "Local", "Nested")
  vector_to_apply=sp_varImpDF_for_graph$model_scale
  library(plyr)
  sp_varImpDF_for_graph$model_scale <- mapvalues(vector_to_apply, from=to_replace, to=replacement)
  
  
  # plot variable importance by model run per species
  EM_vi<-ggplot(sp_varImpDF_for_graph, aes(expl.name, varImp, fill = expl.name)) + 
    geom_col() + facet_wrap(~ model_scale) + ylim(0, 1) +
    scale_fill_manual(values = cb_palette) + 
    labs(x = NULL, y = 'Mean Variable Importance') + ggtitle(paste(sp.name)) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(emVI_fold, sp.name, "_mean_variable_importance_bar_plot.tiff"))
  
  ##########################
  #no nested
  # plot variable importance by model run per species
  sp_varImpDF_for_graph_GL=sp_varImpDF_for_graph[sp_varImpDF_for_graph$model_scale!="Nested",]
  EM_vi<-ggplot(sp_varImpDF_for_graph_GL, aes(expl.name, varImp, fill = expl.name)) + 
    geom_col() + facet_wrap(~ model_scale) + ylim(0, 1) +
    scale_fill_manual(values = cb_palette) + 
    labs(x = NULL, y = 'Mean Variable Importance') + ggtitle(paste(sp.name)) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(emVI_fold, sp.name, "_mean_variable_importance_bar_plot_GL.tiff"))
  
  
}
dput(names(SS_results_DF))
SS_results_DF$weighted_SS_GL_diff=SS_results_DF$SS_GL_diff*SS_results_DF$GL_meanImp
write.csv(SS_results_DF, paste0(rc_fold,"response_deviations.csv"), row.names = F)
#View(SS_results_DF)

SS_results_DF_spp=SS_results_DF[,c("species", "weighted_SS_GL_diff")]
SS_results_DF_spp=aggregate(SS_results_DF_spp[,2], by=list(SS_results_DF_spp$species), FUN=mean)
names(SS_results_DF_spp)=c("species", "weighted_SS_GL_diff")
write.csv(SS_results_DF_spp, paste0(rc_fold,"species_mean_response_deviations.csv"), row.names = F)
#View(SS_results_DF_spp)

##### END MEAN RESPONSE CURVES #####
