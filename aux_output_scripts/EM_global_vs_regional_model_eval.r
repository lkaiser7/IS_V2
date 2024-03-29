#wd="D:/projects/Invasives_modeling/results/xfirst round of results/main_results/"
setwd(rootDir)
dir.create("combined_results/model_eval_metric/", showWarnings = F, recursive = T)
#eval_metrics=c("TSS", "ROC", "KAPPA") #

eval_stat = eval_stats[1]
for (eval_stat in eval_stats){ #global_notHI regional_HI nested_HI

  all_sp_evalDF_short=read.csv(paste0(rootDir, "combined_results/all_EM_eval_short.csv"))
  all_sp_evalDF_short=all_sp_evalDF_short[,c(-1, -3)]
  names(all_sp_evalDF_short)=c("Species", "Global", "Regional", "Nested", "Skill_deviation")
  plot(all_sp_evalDF_short$Global, all_sp_evalDF_short$Regional)

  geom.text.size = 2
  theme.size = (14/5) * geom.text.size
  
  library(ggplot2)
  a=ggplot(all_sp_evalDF_short, aes(x=Global, y=Regional)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=all_sp_evalDF_short$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+geom_smooth(method = "lm", se = TRUE)+xlab("Global model skill")+ylab("Regional model skill")
  a
  cor(all_sp_evalDF_short$Global, all_sp_evalDF_short$Regional)
  
  tiff_name=paste0("combined_results/model_eval_metric/EM_eval_metric_comparison_EM_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  

  ######################
  #now compare match in model skill and variable importance
  all_sp_evalDF_short2=all_sp_evalDF_short[,c("Species", "Global", "Regional", "Skill_deviation")]
  species_var_imp_deviations_df=read.csv(paste0(rootDir, "combined_results/all_sp_varImpDF_short_GL.csv"))
  species_var_imp_deviations_df=species_var_imp_deviations_df[,c("current_spp_name", "SS")]
  names(species_var_imp_deviations_df)=c("Species", "varImp_deviation")
  #View(species_var_imp_deviations_df)
  skill_vs_varImp= merge(all_sp_evalDF_short2, species_var_imp_deviations_df, by="Species")
  #View(skill_vs_varImp)
  #names(skill_vs_varImp)
  cor(skill_vs_varImp$Skill_deviation, skill_vs_varImp$varImp_deviation)
  a=ggplot(skill_vs_varImp, aes(x=varImp_deviation, y=Skill_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_vs_varImp$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Variable importance deviation between global and regional scales")+
    ylab("Skill deviation between global and regional models")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/EM_skill_vs_varImp_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  cor(skill_vs_varImp$Global, skill_vs_varImp$varImp_deviation)
  cor(skill_vs_varImp$Regional, skill_vs_varImp$varImp_deviation)
  a=ggplot(skill_vs_varImp, aes(x=varImp_deviation, y=Global)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_vs_varImp$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Variable importance deviation between global and regional scales")+
    ylab("Global model skill")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/EM_global_skill_vs_varImp_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  a=ggplot(skill_vs_varImp, aes(x=varImp_deviation, y=Regional)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_vs_varImp$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Variable importance deviation between global and regional scales")+
    ylab("Regional model skill")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/EM_regional_skill_vs_varImp_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  #now compare model skill and response curve deviation
  rc_fold<-paste0(rootDir, "combined_results/EM_mean_response_curves/")
  # rc_dev_DF=read.csv(paste0(rc_fold,"scld_species_mean_response_deviations.csv"))
  # rc_dev_DF=rc_dev_DF[,c("current_spp_name", "scld_weighted_SS_GL_diff")]
  # names(rc_dev_DF)=c("Species", "RC_dev")
  rc_dev_DF=read.csv(paste0(rc_fold,"all_response_deviations_indicators.csv"))
  names(rc_dev_DF)=c("Species", "species", "rcDev_scld_meanWgt", "rcDev_scld_noWgt", 
                     "rcDev_noscld_noWgt", "rcDev_scld_maxWgt")
  rc_dev_DF=rc_dev_DF[,-2]
  skill_varImp_RC_DF= merge(skill_vs_varImp, rc_dev_DF, by="Species")
  
  #now merge standard niche overlap
  niche_overlap_DF=read.csv(paste0("combined_results/combined_maps/niche_overlap_metric_DF.csv"))
  niche_overlap_DF=niche_overlap_DF[,c("current_spp_name", "GL_clipSuit_I", "GL_Suit_I", "GL_clipSuit_D", "GL_Suit_D")]
  names(niche_overlap_DF)=c("Species", "GL_clipSuit_I", "GL_Suit_I", "GL_clipSuit_D", "GL_Suit_D")
  skill_varImp_RC_DF= merge(skill_varImp_RC_DF, niche_overlap_DF, by="Species")
  
  #merge species introduction year
  spp_establishment_DF=read.csv(paste0(rootDir, "data/spp_establishment_year.csv"))
  spp_establishment_DF$species=replace_spp_names(spp_establishment_DF$species)
  names(spp_establishment_DF)=c("Species", "establishment")
  skill_varImp_RC_DF= merge(skill_varImp_RC_DF, spp_establishment_DF, by="Species")
  
  range(skill_varImp_RC_DF$GL_Suit_I, na.rm=T)
  
  dput(names(skill_varImp_RC_DF))
  library(plotly)
  library(GGally)
  p <- ggpairs(skill_varImp_RC_DF[,c("rcDev_scld_meanWgt", "rcDev_scld_noWgt", "rcDev_noscld_noWgt", 
                                     "rcDev_scld_maxWgt", "GL_clipSuit_I", "GL_Suit_I", "GL_clipSuit_D", 
                                     "GL_Suit_D", "establishment")], title="") 
  ggplotly(p)
  
  
  #establishment year: max weight better
  p <- ggpairs(skill_varImp_RC_DF[,c("rcDev_scld_meanWgt", "rcDev_scld_noWgt", "rcDev_noscld_noWgt", 
                                     "rcDev_scld_maxWgt", "establishment")], title="") 
  ggplotly(p)
  
  #establishment year: max weight better
  p <- ggpairs(skill_varImp_RC_DF[,c("rcDev_scld_meanWgt", "rcDev_scld_noWgt", "rcDev_noscld_noWgt", 
                                     "rcDev_scld_maxWgt", "GL_Suit_I", "GL_Suit_D")], title="") 
  ggplotly(p)
  
  p <- ggpairs(skill_varImp_RC_DF[,-1], title="") 
  ggplotly(p)
  
  
  library(ggcorrplot)
  ggcorrplot(cor(skill_varImp_RC_DF[, -1]))
  #View(skill_varImp_RC_DF)
  
  plot(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$rcDev_scld_maxWgt)
  cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$rcDev_scld_maxWgt)
  
  plot(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Skill_deviation)
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Skill_deviation)
  
  plot(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Global)
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Global)

  plot(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Regional)
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Regional)
  
  cor(skill_varImp_RC_DF$varImp, skill_varImp_RC_DF$Skill_deviation)
  
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$Skill_deviation)
  a=ggplot(skill_varImp_RC_DF, aes(x=rcDev_scld_maxWgt, y=Skill_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Response deviation between global and regional scales")+
    ylab("Skill deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/EM_skill_vs_RC_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$rcDev_scld_maxWgt)
  a=ggplot(skill_varImp_RC_DF, aes(x=rcDev_scld_maxWgt, y=varImp_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Variable importance deviation between global and regional scales")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/varImpDev_vs_RC_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  # cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$GL_clipSuit_I)
  # a=ggplot(skill_varImp_RC_DF, aes(x=GL_clipSuit_I, y=rcDev_scld_maxWgt)) + 
  #   geom_point(aes(size=1.25)) +
  #   geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
  #   theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
  #   ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  # a
  # tiff_name=paste0("combined_results/model_eval_metric/RC_SSdeviation_vs_warrenIoverlap_", eval_stat, ".tiff")
  # ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")

  
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$GL_Suit_I)
  a=ggplot(skill_varImp_RC_DF, aes(x=GL_Suit_I, y=rcDev_scld_maxWgt)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/RC_SSdeviation_vs_warrenIoverlap_unclipped_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  # cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$GL_Suit_I)
  # a=ggplot(skill_varImp_RC_DF, aes(x=GL_Suit_I, y=varImp_deviation)) + 
  #   geom_point(aes(size=1.25)) +
  #   geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
  #   theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
  #   ylab("Variable importance deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  # a
  # tiff_name=paste0("combined_results/model_eval_metric/varImpdeviation_vs_warrenIoverlap_unclipped_", eval_stat, ".tiff")
  # ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  ####################
  #Shroner's D
  # cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$GL_clipSuit_D)
  # a=ggplot(skill_varImp_RC_DF, aes(x=GL_clipSuit_D, y=rcDev_scld_maxWgt)) + 
  #   geom_point(aes(size=1.25)) +
  #   geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
  #   theme(legend.position="none")+xlab("Niche overlap (Schoener's D) between global and reg. projections")+
  #   ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  # a
  # tiff_name=paste0("combined_results/model_eval_metric/RC_SSdeviation_vs_SchoenerD_overlap_", eval_stat, ".tiff")
  # ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$GL_Suit_D)
  a=ggplot(skill_varImp_RC_DF, aes(x=GL_Suit_D, y=rcDev_scld_maxWgt)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Niche overlap (Schoener's D) between global and reg. projections")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/RC_SSdeviation_vs_SchoenerD_overlap_unclipped_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  # cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$GL_Suit_D)
  # a=ggplot(skill_varImp_RC_DF, aes(x=GL_clipSuit_I, y=varImp_deviation)) + 
  #   geom_point(aes(size=1.25)) +
  #   geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
  #   theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
  #   ylab("Variable importance deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  # a
  # tiff_name=paste0("combined_results/model_eval_metric/varImpdeviation_vs_warrenIoverlap_", eval_stat, ".tiff")
  # ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  ########################
  # cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$GL_Suit_I)
  # a=ggplot(skill_varImp_RC_DF, aes(x=GL_Suit_I, y=varImp_deviation)) + 
  #   geom_point(aes(size=1.25)) +
  #   geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
  #   theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
  #   ylab("Variable importance deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  # a
  # tiff_name=paste0("combined_results/model_eval_metric/varImpdeviation_vs_warrenIoverlap_unclipped_", eval_stat, ".tiff")
  # ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$rcDev_scld_maxWgt, skill_varImp_RC_DF$establishment)
  a=ggplot(skill_varImp_RC_DF, aes(x=establishment, y=rcDev_scld_maxWgt)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Establishment year")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/establishment_vs_RC_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  graphics.off()
  
  # cor(skill_varImp_RC_DF$GL_Suit_I, skill_varImp_RC_DF$establishment)
  # a=ggplot(skill_varImp_RC_DF, aes(x=establishment, y=GL_Suit_I)) + 
  #   geom_point(aes(size=1.25)) +
  #   geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
  #   theme(legend.position="none")+xlab("Establishment year")+
  #   ylab("Niche overlap (Warren's I) between global and reg. projections")+geom_smooth(method = "lm", se = TRUE)
  # a
  # tiff_name=paste0("combined_results/model_eval_metric/establishment_vs_warrensI_", eval_stat, ".tiff")
  # ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$GL_Suit_D, skill_varImp_RC_DF$establishment)
  a=ggplot(skill_varImp_RC_DF, aes(x=establishment, y=GL_Suit_D)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Establishment year")+
    ylab("Niche overlap (Shoener's D) between global and reg. projections")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/establishment_vs_shoenerD_unclipped_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  write.csv(skill_varImp_RC_DF, paste0("combined_results/model_eval_metric/skill_varImp_RC_DF_", eval_stat, ".csv"), row.names = F)
}

