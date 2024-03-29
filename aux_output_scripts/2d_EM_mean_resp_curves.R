### mean ensemble modeling curves ###
##################
##### SET UP #####
##################

### FROM MASTER SCRIPT ###
# set working directory to main analysis folder
setwd(rootDir)

# list model scales
all_mod_scale = c("global_notHI", "regional_HI", "nested_HI")#c("regional", "global", "nested")
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
model_scales=c("global_notHI", "regional_HI", "nested_HI")
model_scale=model_scales[1]
for (model_scale in model_scales){
  tmp_project_path=paste0(rootDir, model_scale, "_models/")
  cat("doing ", model_scale, "\n")
  for (sp_nm in all_sp_nm){
    cat("doing ", sp_nm, "\n")
    current_sp_nm=replace_spp_names(sp_nm)
    
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

all_sp_varImpDF=cbind(all_sp_varImpDF, current_spp_name=replace_spp_names(all_sp_varImpDF$species))

dir.create(paste0(rootDir, "combined_results/"), showWarnings = F)
write.csv(all_sp_varImpDF, paste0(rootDir, "combined_results/all_EM_var_imp.csv"))

#reformat
library(reshape2)
all_sp_varImpDF_short=dcast(all_sp_varImpDF, species ~ pred + scale, value.var="varImp")
all_sp_varImpDF_short_GL=dcast(all_sp_varImpDF[all_sp_varImpDF$scale!="nested_HI",], species ~ pred + scale, value.var="varImp")

all_sp_varImpDF_S_short_GL=dcast(all_sp_varImpDF[all_sp_varImpDF$scale!="nested_HI",], pred + species ~ scale, value.var="varImp")
all_sp_varImpDF_S_short_GL$SS=(all_sp_varImpDF_S_short_GL$global_notHI-all_sp_varImpDF_S_short_GL$regional_HI)^2
all_sp_varImpDF_S_short_GL=all_sp_varImpDF_S_short_GL[,c("pred", "species", "SS")]
all_sp_varImpDF_SS_short_GL=dcast(all_sp_varImpDF_S_short_GL, species ~ pred, value.var="SS")
all_sp_varImpDF_SS_short_GL$SS=apply(all_sp_varImpDF_SS_short_GL[,-1], 1, sum)
all_sp_varImpDF_SS_short_GL=all_sp_varImpDF_SS_short_GL[, c("species", "SS")]
all_sp_varImpDF_short_GL=merge(all_sp_varImpDF_short_GL, all_sp_varImpDF_SS_short_GL, by="species")

all_sp_varImpDF_short=cbind(current_spp_name=replace_spp_names(all_sp_varImpDF_short$species), all_sp_varImpDF_short)
all_sp_varImpDF_short_GL=cbind(current_spp_name=replace_spp_names(all_sp_varImpDF_short_GL$species), all_sp_varImpDF_short_GL)

write.csv(all_sp_varImpDF_short, paste0(rootDir, "combined_results/all_sp_varImpDF_short.csv"))
write.csv(all_sp_varImpDF_short_GL, paste0(rootDir, "combined_results/all_sp_varImpDF_short_GL.csv"))

#####################################################
######################
#second  collect EM model eval
sp.name = all_sp_nm[1]
model_scales=c("global_notHI", "regional_HI", "nested_HI")
model_scale=model_scales[1]
for (model_scale in model_scales){
  tmp_project_path=paste0(rootDir, model_scale, "_models/")
  cat("doing ", model_scale, "\n")
  for (sp_nm in all_sp_nm){
    cat("doing ", sp_nm, "\n")
    
    # replace species naming convention of "_" with "." 
    sp.name_period=str_replace_all(sp_nm, "_", ".")
    sp_dir = paste0(sp.name_period, "/")
    
    FileName00=paste0(tmp_project_path, sp_dir, sp.name_period, "_",eval_stats,"_EM.csv")
    sp_evalDF=read.csv(FileName00)
    sp_evalDF=sp_evalDF[,-1]
    sp_evalDF=as.data.frame(t(matrix(sp_evalDF)))
    names(sp_evalDF)=c("eval", "cutoff", "sensitivity", "specificity")
    sp_evalDF$scale=model_scale
    sp_evalDF$species=sp_nm
    if (sp_nm == all_sp_nm[1] & model_scale==model_scales[1]){
      all_sp_evalDF=sp_evalDF
    }else{
      all_sp_evalDF=rbind(all_sp_evalDF, sp_evalDF)
    }
  }  
}
all_sp_evalDF=cbind(all_sp_evalDF, current_spp_name=replace_spp_names(all_sp_evalDF$species))
#View(all_sp_evalDF)
write.csv(all_sp_evalDF, paste0(rootDir, "combined_results/all_EM_eval.csv"))

library(reshape2)
all_sp_evalDF_short=dcast(all_sp_evalDF[,-2], species ~ scale, value.var="eval")
all_sp_evalDF_short$GL_diff=abs(all_sp_evalDF_short$global_notHI-all_sp_evalDF_short$regional_HI)
all_sp_evalDF_short=cbind(current_spp_name=replace_spp_names(all_sp_evalDF_short$species), all_sp_evalDF_short)

write.csv(all_sp_evalDF_short, paste0(rootDir, "combined_results/all_EM_eval_short.csv"))

all_sp_sensDF_short=dcast(all_sp_evalDF[,-2], species ~ scale, value.var="sensitivity")
all_sp_sensDF_short$GL_diff=abs(all_sp_sensDF_short$global_notHI-all_sp_sensDF_short$regional_HI)
all_sp_sensDF_short=cbind(current_spp_name=replace_spp_names(all_sp_sensDF_short$species), all_sp_sensDF_short)
write.csv(all_sp_sensDF_short, paste0(rootDir, "combined_results/all_EM_sensitivity_short.csv"))

all_sp_specDF_short=dcast(all_sp_evalDF[,-2], species ~ scale, value.var="specificity")
all_sp_specDF_short$GL_diff=abs(all_sp_specDF_short$global_notHI-all_sp_specDF_short$regional_HI)
all_sp_specDF_short=cbind(current_spp_name=replace_spp_names(all_sp_specDF_short$species), all_sp_specDF_short)
write.csv(all_sp_specDF_short, paste0(rootDir, "combined_results/all_EM_specificity_short.csv"))

#combined
tmp1=all_sp_evalDF_short[,c(-1, -6)]
names(tmp1)[2:4]=paste0("eval_", names(tmp1)[2:4])

tmp2=all_sp_sensDF_short[,c(-1, -6)]
names(tmp2)[2:4]=paste0("sens_", names(tmp2)[2:4])

tmp3=all_sp_specDF_short[,c(-1, -6)]
names(tmp3)[2:4]=paste0("spec_", names(tmp3)[2:4])

tmp1=merge(tmp1, tmp2, by="species")
tmp1=merge(tmp1, tmp3, by="species")
tmp1=cbind(current_spp_name=replace_spp_names(tmp1$species), tmp1)

write.csv(tmp1, paste0(rootDir, "combined_results/all_EM_eval_table.csv"), row.names = F)
apply(tmp1[,c(-1, -2)], 2, mean, na.rm=T)

#####################################################
#####################################################
# create folder in output directory for mean response curve plots
rc_fold<-paste0(rootDir, "combined_results/EM_mean_response_curves/")
dir.create(rc_fold, showWarnings = FALSE, recursive = T) 

emVI_fold<-paste0(rootDir, "combined_results/EM_mean_VariImp_plots/")
dir.create(emVI_fold, showWarnings = FALSE, recursive = T) 
dir.create(paste0(rc_fold, "individual_plots/"), showWarnings = FALSE, recursive = T) 
dir.create(paste0(rc_fold, "GL_scld/"), showWarnings = FALSE, recursive = T) 

# loop through each species
s=1
for(s in 1:length(all_sp_nm)){ # set s = 1 for debugging
  # select single species name
  sp.name<-all_sp_nm[s]
  current_sp_nm=replace_spp_names(sp.name)
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
    
    scld_rp.dat2=rp.dat2
    bio = var_names[1]
    for (bio in var_names){
      jnk=scld_rp.dat2[scld_rp.dat2$expl.name==bio,"pred.val"]
      scld_rp.dat2[scld_rp.dat2$expl.name==bio,"pred.val"]=c(scale(jnk))
    }
    # save tables per model scale run
    if(mod_scale == "regional_HI"){ #"global_notHI", "regional_HI", "nested_HI"
      regional_tab<-rp.dat2
      scld_regional_tab<-scld_rp.dat2
    }else if(mod_scale == "global_notHI"){
      global_tab<-rp.dat2
      scld_global_tab<-scld_rp.dat2
    }else{
      nested_tab<-rp.dat2
      scld_nested_tab<-scld_rp.dat2
    }
    
    # then merge to make a the panel figure
    
    # add model scale column
    rp.dat2$model_scale=mod_scale
    scld_rp.dat2$model_scale=mod_scale
    # merge data for plotting
    if(m == 1){
      all.rp.dat2<-rp.dat2
      scld_all.rp.dat2<-scld_rp.dat2
      #all_var_imp<-sp_var_imp
    }else{
      all.rp.dat2<-rbind(all.rp.dat2, rp.dat2)
      scld_all.rp.dat2<-rbind(scld_all.rp.dat2, scld_rp.dat2)
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
  
  scld_plot_data<-merge(scld_all.rp.dat2, sp_varImpDF, by =  c("expl.name", "model_scale"), all = T)
  scld_plot_data$line_size<-scld_plot_data$varImp*3
  
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
  #View(mean_var_imp_df)
  names(mean_var_imp_df)=c("expl.name", "GvarImp", "LvarImp", "NvarImp")
  # mean_var_imp_df$GL_Imp=min(mean_var_imp_df$GvarImp, mean_var_imp_df$LvarImp)
  # mean_var_imp_df$GN_Imp=min(mean_var_imp_df$GvarImp,mean_var_imp_df$NvarImp)
  # mean_var_imp_df$LN_Imp=min(mean_var_imp_df$LvarImp,mean_var_imp_df$NvarImp)
  
  mean_var_imp_df$GL_Imp=(mean_var_imp_df$GvarImp+mean_var_imp_df$LvarImp)/2
  mean_var_imp_df$GN_Imp=(mean_var_imp_df$GvarImp+mean_var_imp_df$NvarImp)/2
  mean_var_imp_df$LN_Imp=(mean_var_imp_df$LvarImp+mean_var_imp_df$NvarImp)/2
  
  mean_var_imp_df$GL_Imp_max=apply(mean_var_imp_df[,c("GvarImp", "LvarImp")],1, max, na.rm=T)
  mean_var_imp_df$GN_Imp_max=apply(mean_var_imp_df[,c("GvarImp", "NvarImp")],1, max, na.rm=T)
  mean_var_imp_df$LN_Imp_max=apply(mean_var_imp_df[,c("LvarImp", "NvarImp")],1, max, na.rm=T)
  
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
    extrapolated_response_DF$GL_diff=(extrapolated_response_DF$global_notHI_models-extrapolated_response_DF$regional_HI_models)^2
    extrapolated_response_DF$GN_diff=(extrapolated_response_DF$global_notHI_models-extrapolated_response_DF$nested_HI_models)^2
    extrapolated_response_DF$LN_diff=(extrapolated_response_DF$regional_HI_models-extrapolated_response_DF$nested_HI_models)^2
    # extrapolated_response_DF$GL_diff=abs(extrapolated_response_DF$global_notHI_models-extrapolated_response_DF$regional_HI_models)
    # extrapolated_response_DF$GN_diff=abs(extrapolated_response_DF$global_notHI_models-extrapolated_response_DF$nested_HI_models)
    # extrapolated_response_DF$LN_diff=abs(extrapolated_response_DF$regional_HI_models-extrapolated_response_DF$nested_HI_models)
    extrapolated_response_DF$scld_global_notHI_models=c(scale(extrapolated_response_DF$global_notHI_models))
    extrapolated_response_DF$scld_regional_HI_models=c(scale(extrapolated_response_DF$regional_HI_models))
    extrapolated_response_DF$scld_nested_HI_models=c(scale(extrapolated_response_DF$nested_HI_models))
    # names(extrapolated_response_DF)=c("expl_vals", "global_notHI_models", "regional_HI_models", "nested_HI_models", 
    #                                   "GL_diff", "GN_diff", "LN_diff", "scld_global_notHI_models", 
    #                                   "scld_regional_HI_models", "scld_nested_HI_models")
    extrapolated_response_DF$scld_GL_diff=(extrapolated_response_DF$scld_global_notHI_models-extrapolated_response_DF$scld_regional_HI_models)^2
    extrapolated_response_DF$scld_GN_diff=(extrapolated_response_DF$scld_global_notHI_models-extrapolated_response_DF$scld_nested_HI_models)^2
    extrapolated_response_DF$scld_LN_diff=(extrapolated_response_DF$scld_regional_HI_models-extrapolated_response_DF$scld_nested_HI_models)^2
    #View(extrapolated_response_DF)
    
    extrapolated_response_DF_tmp=extrapolated_response_DF[,c("expl_vals", "global_notHI_models", "regional_HI_models", "nested_HI_models")]
    extrapolated_response_DF_tmp=melt(extrapolated_response_DF_tmp,id.vars = "expl_vals")
    extrapolated_response_DF_tmp$pred=bio
    extrapolated_response_DF_tmp$species=sp.name
    ggplot(extrapolated_response_DF_tmp, aes(x = expl_vals, y = value, colour = variable))+geom_line()
    #View(extrapolated_response_DF_tmp)
    
    ###create_extrapolated scaled response df
    scld_extrapolated_response_DF_tmp=extrapolated_response_DF[,c("expl_vals", "scld_global_notHI_models", "scld_regional_HI_models", "scld_GL_diff")]
    names(scld_extrapolated_response_DF_tmp)=c("expl_vals", "global_notHI_models", "regional_HI_models", "squared_distance")
    scld_extrapolated_response_DF_tmp=melt(scld_extrapolated_response_DF_tmp,id.vars = "expl_vals")
    scld_extrapolated_response_DF_tmp$pred=bio
    scld_extrapolated_response_DF_tmp$species=sp.name
    a=ggplot(scld_extrapolated_response_DF_tmp, aes(x = expl_vals, y = value, colour = variable))+geom_line()+xlab(bio)+ylab("Scaled suitability")
    ggsave(a, filename = paste0(rc_fold, "individual_plots/", sp.name, "_", bio, "_EM_scaled_GL_resp_curve.tiff"), compression = "lzw")
    
    # extrapolated_response_DF_tmp=extrapolated_response_DF[,c("expl_vals", "GL_diff", "GN_diff", "LN_diff")]
    # extrapolated_response_DF_tmp=melt(extrapolated_response_DF_tmp,id.vars = "expl_vals")
    # ggplot(extrapolated_response_DF_tmp, aes(x = expl_vals, y = value, colour = variable))+geom_line()
    
    # scld_extrapolated_response_DF_tmp=extrapolated_response_DF[,c("expl_vals", "scld_GL_diff", "scld_GN_diff", "scld_LN_diff")]
    # scld_extrapolated_response_DF_tmp=melt(scld_extrapolated_response_DF_tmp,id.vars = "expl_vals")
    # ggplot(scld_extrapolated_response_DF_tmp, aes(x = expl_vals, y = value, colour = variable))+geom_line()
    
    # SS_G_L_diff=sum(extrapolated_response_DF$G_L_diff^2)
    # SS_G_N_diff=sum(extrapolated_response_DF$G_N_diff^2)
    # SS_L_N_diff=sum(extrapolated_response_DF$L_N_diff^2)
    
    SS_GL_diff=sum(extrapolated_response_DF$GL_diff)
    SS_GN_diff=sum(extrapolated_response_DF$GN_diff)
    SS_LN_diff=sum(extrapolated_response_DF$LN_diff)
    
    scld_SS_GL_diff=sum(extrapolated_response_DF$scld_GL_diff)
    scld_SS_GN_diff=sum(extrapolated_response_DF$scld_GN_diff)
    scld_SS_LN_diff=sum(extrapolated_response_DF$scld_LN_diff)
    
    bio_mean_var_imp_df=mean_var_imp_df[mean_var_imp_df$expl.name==bio,-1]
    results_row=data.frame(model=mod_type, species=sp.name, bio, SS_GL_diff, SS_GN_diff, SS_LN_diff, scld_SS_GL_diff, scld_SS_GN_diff, scld_SS_LN_diff) 
    results_row=cbind(results_row, bio_mean_var_imp_df)
    
    if (s==1 & bio==var_names[1]){
      SS_results_DF=results_row
      all_extrapolated_response_DF=extrapolated_response_DF_tmp
      all_scld_extrapolated_response_DF=scld_extrapolated_response_DF_tmp
    }else{
      SS_results_DF=rbind(SS_results_DF, results_row)
      all_extrapolated_response_DF=rbind(all_extrapolated_response_DF, extrapolated_response_DF_tmp)
      all_scld_extrapolated_response_DF=rbind(all_scld_extrapolated_response_DF, scld_extrapolated_response_DF_tmp)
    }
  }

  #View(SS_results_DF)
  # ggplot(plot_data,
  #        aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
  #   geom_line(lwd = plot_data$line_size) + #geom_line(lwd = plot_data$pred.mean) + 
  #   facet_wrap(~ expl.name, scales = 'free_x') +
  #   # coord_cartesian(xlim = c(0, 150)) +  # sets manual x scale limits
  #   facet_wrap_custom(~ expl.name, scales = "free", scale_overrides = list(
  #     scale_override(1, scale_x_continuous(limits = c(0, 50))),
  #     scale_override(3, scale_x_continuous(limits = c(0, 150))),
  #     scale_override(4, scale_x_continuous(limits = c(0, 40)))
  #   )) +
  #   labs(x = '', y = 'probability of occurrence', colour = 'model scale') + 
  #   # scale_color_brewer(type = 'qual', palette = 4) + 
  #   scale_color_manual(values = cb_palette) + 
  #   ggtitle(paste(sub("_", " ", sp.name))) + 
  #   theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  # 
  # ggsave(filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve.tiff"))
  #View(scld_plot_data)
  
  # ggplot(scld_plot_data,
  #        aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
  #   geom_line(lwd = plot_data$line_size) + #geom_line(lwd = plot_data$pred.mean) + 
  #   facet_wrap(~ expl.name, scales = 'free_x') +
  #   # coord_cartesian(xlim = c(0, 150)) +  # sets manual x scale limits
  #   # facet_wrap_custom(~ expl.name, scales = "free", scale_overrides = list(
  #   #   scale_override(1, scale_x_continuous(limits = c(0, 50))),
  #   #   scale_override(3, scale_x_continuous(limits = c(0, 150))),
  #   #   scale_override(4, scale_x_continuous(limits = c(0, 40)))
  #   # )) +
  #   labs(x = '', y = 'Scaled response', colour = 'model scale') + 
  #   # scale_color_brewer(type = 'qual', palette = 4) + 
  #   scale_color_manual(values = cb_palette) + 
  #   ggtitle(current_sp_nm) + 
  #   theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  # 
  # ggsave(filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve_scaled.tiff"), compression = "lzw")
  
  ########
  library(ggpubr)
  leg_pos='none'
  bio = var_names[1]
  plot_list=list()
  ggp=1
  for (bio in var_names){
    tmp_scld_plot_data=scld_plot_data[scld_plot_data$expl.name==bio,]
    tmp_plot_data=plot_data[scld_plot_data$expl.name==bio,]
    if (bio==var_names[length(var_names)]){
      leg_pos='none' 
    }
    tmpplot=ggplot(tmp_scld_plot_data,
           aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
      geom_line(lwd = tmp_plot_data$line_size) +
      labs(x = bio, y = '', colour = 'model scale') + 
      scale_color_manual(values = cb_palette, name=NULL) + 
      # ggtitle(bio) + 
      theme_minimal() + theme(legend.position = leg_pos, plot.title = element_text(hjust = 0.5))+ylab("Scaled suitability")
    plot_list[[ggp]]=tmpplot
    #assign(paste0(bio, "_plot"), tmpplot)
    ggp=ggp+1
  }
  a=ggarrange(plotlist=plot_list)
  ggsave(plot = a, filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve_scaled.tiff"), compression = "lzw", height=7, width=7, units = "in", dpi = 300)
  
  ##############################
  #no nested:
  #View(plot_data)
  
  # plot_data_no_nested=plot_data[plot_data$model_scale!="nested_HI_models",]
  # 
  # ggplot(plot_data_no_nested,
  #        aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
  #   geom_line(lwd = plot_data_no_nested$line_size) + #geom_line(lwd = plot_data$pred.mean) + 
  #   facet_wrap(~ expl.name, scales = 'free_x') +
  #   # coord_cartesian(xlim = c(0, 150)) +  # sets manual x scale limits
  #   facet_wrap_custom(~ expl.name, scales = "free", scale_overrides = list(
  #     scale_override(1, scale_x_continuous(limits = c(0, 50))),
  #     scale_override(3, scale_x_continuous(limits = c(0, 150))),
  #     scale_override(4, scale_x_continuous(limits = c(0, 40)))
  #   )) +
  #   labs(x = '', y = 'probability of occurrence', colour = 'model scale') + 
  #   # scale_color_brewer(type = 'qual', palette = 4) + 
  #   scale_color_manual(values = cb_palette) + 
  #   ggtitle(current_sp_nm) + 
  #   theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  # 
  # ggsave(filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve_GL.tiff"), compression = "lzw")
  # 
  ###############
  plot_data_no_nested=plot_data[plot_data$model_scale!="nested_HI_models",]
  library(ggpubr)
  leg_pos='none'
  bio = var_names[1]
  plot_list=list()
  ggp=1
  for (bio in var_names){
    tmp_plot_data=plot_data_no_nested[plot_data_no_nested$expl.name==bio,]
    if (bio==var_names[length(var_names)]){
      leg_pos='none' 
    }
    tmpplot=ggplot(tmp_plot_data,
                   aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
      geom_line(lwd = tmp_plot_data$line_size) +
      labs(x = bio, y = '', colour = 'model scale') + 
      scale_color_manual(values = cb_palette, name=NULL) + 
      # ggtitle(bio) + 
      theme_minimal() + theme(legend.position = leg_pos, plot.title = element_text(hjust = 0.5))+ylab("Scaled suitability")
    plot_list[[ggp]]=tmpplot
    #assign(paste0(bio, "_plot"), tmpplot)
    ggp=ggp+1
  }
  a=ggarrange(plotlist=plot_list)
  ggsave(plot = a, filename = paste0(rc_fold, sp.name, "_EM_mean_resp_curve_GL.tiff"), compression = "lzw", height=7, width=7, units = "in", dpi = 300)
  
  #######
  #scaled
  # scld_plot_data_no_nested=scld_plot_data[scld_plot_data$model_scale!="nested_HI_models",]
  # ggplot(scld_plot_data_no_nested,
  #        aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
  #   geom_line(lwd = plot_data_no_nested$line_size) + #geom_line(lwd = plot_data$pred.mean) + 
  #   facet_wrap(~ expl.name, scales = 'free_x') +
  #   # coord_cartesian(xlim = c(0, 150)) +  # sets manual x scale limits
  #   facet_wrap_custom(~ expl.name, scales = "free", scale_overrides = list(
  #     scale_override(1, scale_x_continuous(limits = c(0, 50))),
  #     scale_override(3, scale_x_continuous(limits = c(0, 150))),
  #     scale_override(4, scale_x_continuous(limits = c(0, 40)))
  #   )) +
  #   labs(x = '', y = 'probability of occurrence', colour = 'model scale') + 
  #   # scale_color_brewer(type = 'qual', palette = 4) + 
  #   scale_color_manual(values = cb_palette) + 
  #   ggtitle(current_sp_nm) + 
  #   theme_minimal() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5))
  # 
  # ggsave(filename = paste0(rc_fold, "GL_scld/", sp.name, "_EM_mean_resp_curve_GL_scld.tiff"), compression = "lzw")
  
  #########
  scld_plot_data_no_nested=scld_plot_data[scld_plot_data$model_scale!="nested_HI_models",]
  library(ggpubr)
  leg_pos='none'
  bio = var_names[1]
  plot_list=list()
  ggp=1
  for (bio in var_names){
    tmp_plot_data=scld_plot_data_no_nested[scld_plot_data_no_nested$expl.name==bio,]
    if (bio==var_names[length(var_names)]){
      leg_pos='none' 
    }
    tmpplot=ggplot(tmp_plot_data,
                   aes(x = expl.val, y = pred.val, colour = model_scale))+#, group = expl.name) + 
      geom_line(lwd = tmp_plot_data$line_size) +
      labs(x = bio, y = '', colour = 'model scale') + 
      scale_color_manual(values = cb_palette, name=NULL) + 
      # ggtitle(bio) + 
      theme_minimal() + theme(legend.position = leg_pos, plot.title = element_text(hjust = 0.5))+ylab("Scaled suitability")
    plot_list[[ggp]]=tmpplot
    #assign(paste0(bio, "_plot"), tmpplot)
    ggp=ggp+1
  }
  a=ggarrange(plotlist=plot_list)
  ggsave(plot = a, filename = paste0(rc_fold, "GL_scld/", sp.name, "_EM_mean_resp_curve_GL_scld.tiff"), compression = "lzw", height=7, width=7, units = "in", dpi = 300)
  
  ########
  #revised var imp plots
  # color blind palette
  cb_palette<-c("#009E73", "#E69F00", "#CC79a7", "#56B4E9")
  sp_varImpDF_for_graph=sp_varImpDF
  
  to_replace=c("global_notHI_models", "regional_HI_models", "nested_HI_models")
  replacement=c("Global", "Regional", "Nested")
  vector_to_apply=sp_varImpDF_for_graph$model_scale
  library(plyr)
  sp_varImpDF_for_graph$model_scale <- mapvalues(vector_to_apply, from=to_replace, to=replacement)
  
  
  # plot variable importance by model run per species
  EM_vi<-ggplot(sp_varImpDF_for_graph, aes(expl.name, varImp, fill = expl.name)) + 
    geom_col() + facet_wrap(~ model_scale) + ylim(0, 1) +
    scale_fill_manual(values = cb_palette) + 
    labs(x = NULL, y = 'Mean Variable Importance') + ggtitle(current_sp_nm) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(emVI_fold, sp.name, "_mean_variable_importance_bar_plot.tiff"), compression = "lzw")
  
  ##########################
  #no nested
  # plot variable importance by model run per species
  sp_varImpDF_for_graph_GL=sp_varImpDF_for_graph[sp_varImpDF_for_graph$model_scale!="Nested",]
  EM_vi<-ggplot(sp_varImpDF_for_graph_GL, aes(expl.name, varImp, fill = expl.name)) + 
    geom_col() + facet_wrap(~ model_scale) + ylim(0, 1) +
    scale_fill_manual(values = cb_palette) + 
    labs(x = NULL, y = 'Mean Variable Importance') + ggtitle(current_sp_nm) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(emVI_fold, sp.name, "_mean_variable_importance_bar_plot_GL.tiff"), compression = "lzw")
  
  
}

#View(SS_results_DF_scld_SS_GL_diff_short)
#################
#DEVIATIONS BASED ON MEAN VAR IMP
SS_results_DF$weighted_SS_GL_diff=SS_results_DF$SS_GL_diff*SS_results_DF$GL_Imp
SS_results_DF$scld_weighted_SS_GL_diff=SS_results_DF$scld_SS_GL_diff*SS_results_DF$GL_Imp
SS_results_DF=cbind(current_spp_name=replace_spp_names(SS_results_DF$species), SS_results_DF)
write.csv(SS_results_DF, paste0(rc_fold,"response_deviations.csv"), row.names = F)
#View(SS_results_DF)

SS_results_DF_spp=SS_results_DF[,c("species", "weighted_SS_GL_diff")]
SS_results_DF_spp=aggregate(SS_results_DF_spp[,2], by=list(SS_results_DF_spp$species), FUN=mean)
names(SS_results_DF_spp)=c("species", "weighted_SS_GL_diff")
SS_results_DF_spp=cbind(current_spp_name=replace_spp_names(SS_results_DF_spp$species), SS_results_DF_spp)
write.csv(SS_results_DF_spp, paste0(rc_fold,"species_mean_response_deviations.csv"), row.names = F)
#View(SS_results_DF_spp)

#############
#scaled and weighted table
SS_results_DF_scld_weighted_SS_GL_diff=SS_results_DF[,c("species", "scld_weighted_SS_GL_diff", "bio")]
library(reshape2)
#View(SS_results_DF_scld_weighted_SS_GL_diff)
SS_results_DF_scld_weighted_SS_GL_diff_short=dcast(SS_results_DF_scld_weighted_SS_GL_diff, species ~ bio, value.var="scld_weighted_SS_GL_diff")

scld_weighted_SS_results_DF_spp=SS_results_DF[,c("species", "scld_weighted_SS_GL_diff")]
scld_weighted_SS_results_DF_spp=aggregate(scld_weighted_SS_results_DF_spp[,2], by=list(scld_weighted_SS_results_DF_spp$species), FUN=mean)
names(scld_weighted_SS_results_DF_spp)=c("species", "scld_weighted_SS_GL_diff")
scld_weighted_SS_results_DF_spp=merge(SS_results_DF_scld_weighted_SS_GL_diff_short, scld_weighted_SS_results_DF_spp, by= "species")
scld_weighted_SS_results_DF_spp=cbind(current_spp_name=replace_spp_names(scld_weighted_SS_results_DF_spp$species), scld_weighted_SS_results_DF_spp)
write.csv(scld_weighted_SS_results_DF_spp, paste0(rc_fold,"scld_weighted_species_mean_response_deviations.csv"), row.names = F)
#View(scld_SS_results_DF_spp)

#############
#scaled and NOT weighted table
SS_results_DF_scld_SS_GL_diff=SS_results_DF[,c("species", "scld_SS_GL_diff", "bio")]
library(reshape2)
#View(SS_results_DF_scld_SS_GL_diff)
SS_results_DF_scld_SS_GL_diff_short=dcast(SS_results_DF_scld_SS_GL_diff, species ~ bio, value.var="scld_SS_GL_diff")

scld_SS_results_DF_spp=SS_results_DF[,c("species", "scld_SS_GL_diff")]
scld_SS_results_DF_spp=aggregate(scld_SS_results_DF_spp[,2], by=list(scld_SS_results_DF_spp$species), FUN=mean)
names(scld_SS_results_DF_spp)=c("species", "scld_SS_GL_diff")
scld_SS_results_DF_spp=merge(SS_results_DF_scld_SS_GL_diff_short, scld_SS_results_DF_spp, by= "species")
scld_SS_results_DF_spp=cbind(current_spp_name=replace_spp_names(scld_SS_results_DF_spp$species), scld_SS_results_DF_spp)
write.csv(scld_SS_results_DF_spp, paste0(rc_fold,"scld_NOT_weighted_species_mean_response_deviations.csv"), row.names = F)
#View(scld_SS_results_DF_spp)

#############
#NOT scaled and NOT weighted table
SS_results_DF_SS_GL_diff=SS_results_DF[,c("species", "SS_GL_diff", "bio")]
library(reshape2)
#View(SS_results_DF_SS_GL_diff)
SS_results_DF_SS_GL_diff_short=dcast(SS_results_DF_SS_GL_diff, species ~ bio, value.var="SS_GL_diff")

SS_results_DF_spp=SS_results_DF[,c("species", "SS_GL_diff")]
SS_results_DF_spp=aggregate(SS_results_DF_spp[,2], by=list(SS_results_DF_spp$species), FUN=mean)
names(SS_results_DF_spp)=c("species", "SS_GL_diff")
SS_results_DF_spp=merge(SS_results_DF_SS_GL_diff_short, SS_results_DF_spp, by= "species")
SS_results_DF_spp=cbind(current_spp_name=replace_spp_names(SS_results_DF_spp$species), SS_results_DF_spp)
write.csv(SS_results_DF_spp, paste0(rc_fold,"NOT_scld_NOT_weighted_species_mean_response_deviations.csv"), row.names = F)
#View(SS_results_DF_spp)

#################
#################
#DEVIATIONS BASED ON MAX VAR IMP
SS_results_DF$weighted_SS_GL_diff=SS_results_DF$SS_GL_diff*SS_results_DF$GL_Imp_max
SS_results_DF$scld_weighted_SS_GL_diff=SS_results_DF$scld_SS_GL_diff*SS_results_DF$GL_Imp_max
SS_results_DF=cbind(current_spp_name=replace_spp_names(SS_results_DF$species), SS_results_DF)
write.csv(SS_results_DF, paste0(rc_fold,"response_deviations_max_imp.csv"), row.names = F)
#View(SS_results_DF)

SS_results_DF_spp=SS_results_DF[,c("species", "weighted_SS_GL_diff")]
SS_results_DF_spp=aggregate(SS_results_DF_spp[,2], by=list(SS_results_DF_spp$species), FUN=mean)
names(SS_results_DF_spp)=c("species", "weighted_SS_GL_diff")
SS_results_DF_spp=cbind(current_spp_name=replace_spp_names(SS_results_DF_spp$species), SS_results_DF_spp)
write.csv(SS_results_DF_spp, paste0(rc_fold,"species_mean_response_deviations_max_imp.csv"), row.names = F)
#View(SS_results_DF_spp)

#############
#scaled and weighted table
SS_results_DF_scld_weighted_SS_GL_diff=SS_results_DF[,c("species", "scld_weighted_SS_GL_diff", "bio")]
library(reshape2)
#View(SS_results_DF_scld_weighted_SS_GL_diff)
SS_results_DF_scld_weighted_SS_GL_diff_short=dcast(SS_results_DF_scld_weighted_SS_GL_diff, species ~ bio, value.var="scld_weighted_SS_GL_diff")

scld_weighted_SS_results_DF_spp=SS_results_DF[,c("species", "scld_weighted_SS_GL_diff")]
scld_weighted_SS_results_DF_spp=aggregate(scld_weighted_SS_results_DF_spp[,2], by=list(scld_weighted_SS_results_DF_spp$species), FUN=mean)
names(scld_weighted_SS_results_DF_spp)=c("species", "scld_weighted_SS_GL_diff")
scld_weighted_SS_results_DF_spp=merge(SS_results_DF_scld_weighted_SS_GL_diff_short, scld_weighted_SS_results_DF_spp, by= "species")
scld_weighted_SS_results_DF_spp=cbind(current_spp_name=replace_spp_names(scld_weighted_SS_results_DF_spp$species), scld_weighted_SS_results_DF_spp)
write.csv(scld_weighted_SS_results_DF_spp, paste0(rc_fold,"scld_weighted_species_mean_response_deviations_max_imp.csv"), row.names = F)
#View(scld_SS_results_DF_spp)

##### END MEAN RESPONSE CURVES #####
scld_weighted_SS_results_DF_spp=read.csv(paste0(rc_fold,"scld_weighted_species_mean_response_deviations.csv"))
scld_NOTweighted_SS_results_DF_spp=read.csv(paste0(rc_fold,"scld_NOT_weighted_species_mean_response_deviations.csv"))
NOTscld_NOTweighted_SS_results_DF_spp=read.csv(paste0(rc_fold,"NOT_scld_NOT_weighted_species_mean_response_deviations.csv"))
scld_MAX_weighted_SS_results_DF_spp=read.csv(paste0(rc_fold,"scld_weighted_species_mean_response_deviations_max_imp.csv"))
View(scld_weighted_SS_results_DF_spp)
View(scld_NOTweighted_SS_results_DF_spp)
View(NOTscld_NOTweighted_SS_results_DF_spp)
View(scld_MAX_weighted_SS_results_DF_spp)
all_response_vars=scld_weighted_SS_results_DF_spp[,c("current_spp_name", "species", "scld_weighted_SS_GL_diff")]
all_response_vars=cbind(all_response_vars, scld_NOTweighted_SS_results_DF_spp$scld_SS_GL_diff, NOTscld_NOTweighted_SS_results_DF_spp$SS_GL_diff, scld_MAX_weighted_SS_results_DF_spp$scld_weighted_SS_GL_diff)
names(all_response_vars)=c("current_spp_name", "species", "scld_meanWgt", "scld_noWgt", "noscld_noWgt", "scld_maxWgt")
#View(all_response_vars)
write.csv(all_response_vars, paste0(rc_fold,"all_response_deviations_indicators.csv"), row.names = F)
pairs(all_response_vars[,3:6], pch = 19)
library(plotly)
library(GGally)
p <- ggpairs(all_response_vars[,3:6], title="correlogram with ggpairs()") 
ggplotly(p)