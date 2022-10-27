#wd="D:/projects/Invasives_modeling/results/xfirst round of results/main_results/"
setwd(rootDir)
dir.create("combined_results/model_eval_metric/", showWarnings = F, recursive = T)
#eval_metrics=c("TSS", "ROC", "KAPPA") #

eval_stat = eval_stats[1]
for (eval_stat in eval_stats){ #global_notHI local_HI nested_HI

  all_sp_evalDF_short=read.csv(paste0(rootDir, "combined_results/all_EM_eval_short.csv"))
  all_sp_evalDF_short=all_sp_evalDF_short[,-1]
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
  species_var_imp_deviations_df=species_var_imp_deviations_df[,c("species", "SS")]
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
  rc_dev_DF=read.csv(paste0(rc_fold,"scld_species_mean_response_deviations.csv"))
  rc_dev_DF=rc_dev_DF[,c("species", "scld_weighted_SS_GL_diff")]
  names(rc_dev_DF)=c("Species", "RC_dev")
  
  skill_varImp_RC_DF= merge(skill_vs_varImp, rc_dev_DF, by="Species")
  
  #now merge standard niche overlap
  niche_overlap_DF=read.csv(paste0("combined_results/combined_maps/niche_overlap_metric_DF.csv"))
  niche_overlap_DF=niche_overlap_DF[,c("sp_nm", "GL_clipSuit_I")]
  names(niche_overlap_DF)=c("Species", "GL_clipSuit_I")
  skill_varImp_RC_DF= merge(skill_varImp_RC_DF, niche_overlap_DF, by="Species")
  
  #merge species introduction year
  spp_establishment_DF=read.csv(paste0(rootDir, "data/spp_establishment_year.csv"))
  names(spp_establishment_DF)=c("Species", "establishment")
  
  skill_varImp_RC_DF= merge(skill_varImp_RC_DF, spp_establishment_DF, by="Species")
  
  library(ggcorrplot)
  ggcorrplot(cor(skill_varImp_RC_DF[, -1]))
  #View(skill_varImp_RC_DF)
  
  plot(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$RC_dev)
  cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$RC_dev)
  
  plot(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Skill_deviation)
  cor(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Skill_deviation)
  
  plot(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Global)
  cor(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Global)

  plot(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Regional)
  cor(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Regional)
  
  cor(skill_varImp_RC_DF$varImp, skill_varImp_RC_DF$Skill_deviation)
  
  cor(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$Skill_deviation)
  a=ggplot(skill_varImp_RC_DF, aes(x=RC_dev, y=Skill_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Response deviation between global and regional scales")+
    ylab("Skill deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/EM_skill_vs_RC_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$RC_dev)
  a=ggplot(skill_varImp_RC_DF, aes(x=RC_dev, y=varImp_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Variable importance deviation between global and regional scales")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/varImpDev_vs_RC_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$GL_clipSuit_I)
  a=ggplot(skill_varImp_RC_DF, aes(x=GL_clipSuit_I, y=RC_dev)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/RC_SSdeviation_vs_warrenIoverlap_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")

  cor(skill_varImp_RC_DF$varImp_deviation, skill_varImp_RC_DF$GL_clipSuit_I)
  a=ggplot(skill_varImp_RC_DF, aes(x=GL_clipSuit_I, y=varImp_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Niche overlap (Warren's I) between global and reg. projections")+
    ylab("Variable importance deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/varImpdeviation_vs_warrenIoverlap_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  cor(skill_varImp_RC_DF$RC_dev, skill_varImp_RC_DF$establishment)
  a=ggplot(skill_varImp_RC_DF, aes(x=establishment, y=RC_dev)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_varImp_RC_DF$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Establishment year")+
    ylab("Response deviation between global and regional scales")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/establishment_vs_RC_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 5, units = "in", compress="lzw")
  
  write.csv(skill_varImp_RC_DF, paste0("combined_results/model_eval_metric/skill_varImp_RC_DF_", eval_stat, ".csv"), row.names = F)
}

# #now compare metrics
# if (length(eval_stats)==3){
#   for (eval_stat in eval_stats){
#     file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_maxent_", eval_stat, ".csv")
#     maxent_file=read.csv(file_name)
#     maxent_file$metric=eval_stat
#     assign(paste0("maxent_", eval_stat), maxent_file)
#     file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_GBM_", eval_stat, ".csv")
#     gbm_file=read.csv(file_name)
#     gbm_file$metric=eval_stat
#     assign(paste0("gbm_", eval_stat), gbm_file)
#   }
#   all_maxent_metrics=rbind(maxent_TSS, maxent_ROC, maxent_KAPPA)
#   all_gbm_metrics=rbind(gbm_TSS, gbm_ROC, gbm_KAPPA)
#   
#   all_maxent_metrics$Species=gsub(pattern="_", replacement=" ", all_maxent_metrics$Species)
#   all_gbm_metrics$Species=gsub(pattern="_", replacement=" ", all_gbm_metrics$Species)
#   
#   a=ggplot(data=all_maxent_metrics, aes(x=Species, y=Regional, fill=metric)) +
#     geom_bar(stat="identity", position=position_dodge()) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.title=element_blank()) + ylab("Regional maxent eval. metrics") +xlab("")
#   a
#   tiff_name=paste0("combined_results/model_eval_metric/Regional_HI_maxent_eval_metric_comparison.tiff")
#   ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
#   
#   a=ggplot(data=all_maxent_metrics, aes(x=Species, y=Global, fill=metric)) +
#     geom_bar(stat="identity", position=position_dodge()) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.title=element_blank()) + ylab("Global maxent eval. metrics") +xlab("")
#   a
#   tiff_name=paste0("combined_results/model_eval_metric/global_notHI_maxent_eval_metric_comparison.tiff")
#   ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
#   
#   a=ggplot(data=all_gbm_metrics, aes(x=Species, y=Regional, fill=metric)) +
#     geom_bar(stat="identity", position=position_dodge()) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.title=element_blank()) + ylab("Regional gbm eval. metrics") +xlab("")
#   a
#   tiff_name=paste0("combined_results/model_eval_metric/Regional_HI_gbm_eval_metric_comparison.tiff")
#   ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
#   
#   a=ggplot(data=all_gbm_metrics, aes(x=Species, y=Global, fill=metric)) +
#     geom_bar(stat="identity", position=position_dodge()) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     theme(legend.title=element_blank()) + ylab("Global gbm eval. metrics") +xlab("")
#   a
#   tiff_name=paste0("combined_results/model_eval_metric/global_notHI_gbm_eval_metric_comparison.tiff")
#   ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
# }
