#wd="D:/projects/Invasives_modeling/results/xfirst round of results/main_results/"
setwd(rootDir)
dir.create("combined_results/model_eval_metric/", showWarnings = F, recursive = T)
#eval_metrics=c("TSS", "ROC", "KAPPA") #

eval_stat = eval_stats[1]
for (eval_stat in eval_stats){ #global_notHI regional_HI nested_HI
  global_notHI_eval_df=read.csv(paste0("global_notHI_models/outputs/all_eval_mat_",eval_stat,".csv"))
  regional_HI_eval_df=read.csv(paste0("regional_HI_models/outputs/all_eval_mat_",eval_stat,".csv"))
  
  calculate_eval_mean_eval_metric=function(metric_df, model_eval_metric="MAXENT.Phillips"){
    metric_df_maxent_mean=metric_df[metric_df$rownames.Spp_eval.==model_eval_metric, ]
    metric_df_maxent_mean=metric_df_maxent_mean[,-2]
    jnk=apply(metric_df_maxent_mean[,-1], 1, mean, na.rm=T)
    mean_eval_metric=data.frame(species=metric_df_maxent_mean[, 1], eval_stat=jnk)
    return(mean_eval_metric)
  }
  
  mean_global_notHI_maxent_eval=calculate_eval_mean_eval_metric(global_notHI_eval_df)
  mean_regional_HI_maxent_eval=calculate_eval_mean_eval_metric(regional_HI_eval_df)

  mean_global_notHI_GBM_eval=calculate_eval_mean_eval_metric(global_notHI_eval_df, model_eval_metric="GBM")
  mean_regional_HI_GBM_eval=calculate_eval_mean_eval_metric(regional_HI_eval_df, model_eval_metric="GBM")

  mean_maxent_eval_df=merge(mean_global_notHI_maxent_eval, mean_regional_HI_maxent_eval, by="species")
  mean_GBM_eval_df=merge(mean_global_notHI_GBM_eval, mean_regional_HI_GBM_eval, by="species")
  
  names(mean_maxent_eval_df)=c("Species", "Global", "Regional")
  names(mean_GBM_eval_df)=c("Species", "Global", "Regional")
  
  mean_allmodels_eval_df=merge(mean_maxent_eval_df, mean_GBM_eval_df, by="Species")
  mean_allmodels_eval_df$Global=apply(mean_allmodels_eval_df[,c("Global.x", "Global.y")], 1, FUN=mean, na.rm=T)
  mean_allmodels_eval_df$Regional=apply(mean_allmodels_eval_df[,c("Regional.x", "Regional.y")], 1, FUN=mean, na.rm=T)
  mean_allmodels_eval_df=mean_allmodels_eval_df[,c("Species", "Global", "Regional")]

    #View(mean_maxent_eval_df)
  file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_maxent_", eval_stat, ".csv")
  write.csv(mean_maxent_eval_df, file_name, row.names = F)
  file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_GBM_", eval_stat, ".csv")
  write.csv(mean_GBM_eval_df, file_name, row.names = F)
  file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_allModels_", eval_stat, ".csv")
  write.csv(mean_allmodels_eval_df, file_name, row.names = F)
  
  plot(mean_maxent_eval_df$Global, mean_maxent_eval_df$Regional)
  plot(mean_GBM_eval_df$Global, mean_GBM_eval_df$Regional)
  
  geom.text.size = 2
  theme.size = (14/5) * geom.text.size
  
  library(ggplot2)
  a=ggplot(mean_maxent_eval_df, aes(x=Global, y=Regional)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=mean_maxent_eval_df$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+geom_smooth(method = "lm", se = TRUE)+xlab("Global model skill")+ylab("Regional model skill")
  a
  tiff_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_maxent_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  library(ggplot2)
  a=ggplot(mean_GBM_eval_df, aes(x=Global, y=Regional)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=mean_GBM_eval_df$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+geom_smooth(method = "lm", se = TRUE)+xlab("Global model skill")+ylab("Regional model skill")
  a
  tiff_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_GBM_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")

  a=ggplot(mean_allmodels_eval_df, aes(x=Global, y=Regional)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=mean_maxent_eval_df$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+geom_smooth(method = "lm", se = TRUE)+xlab("Global model skill")+ylab("Regional model skill")
  a
  tiff_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_allModels_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  ######################
  #now compare match in model skill and variable importance
  #mean_allmodels_eval_df$skill_deviation=abs(mean_allmodels_eval_df$Global-mean_allmodels_eval_df$Regional)
  mean_allmodels_eval_df$skill_deviation=(mean_allmodels_eval_df$Global-mean_allmodels_eval_df$Regional)^2
  species_var_imp_deviations_df=read.csv("combined_results/mean_VariImp_plots/mean_deviation_in_global_vs_regional_variable_importance.csv")
  #View(species_var_imp_deviations_df)
  skill_vs_varImp= merge(mean_allmodels_eval_df, species_var_imp_deviations_df, by.x="Species", by.y="species")
  #View(skill_vs_varImp)
  #names(skill_vs_varImp)
  cor(skill_vs_varImp$skill_deviation, skill_vs_varImp$varImp_deviation)
  a=ggplot(skill_vs_varImp, aes(x=varImp_deviation, y=skill_deviation)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=skill_vs_varImp$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")+xlab("Deviation between global and regional model variable importance")+
    ylab("Deviation between global and regional model skill")+geom_smooth(method = "lm", se = TRUE)
  a
  tiff_name=paste0("combined_results/model_eval_metric/skill_vs_varImp_SSdeviation_", eval_stat, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
}

#now compare metrics
if (length(eval_stats)==3){
  for (eval_stat in eval_stats){
    file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_maxent_", eval_stat, ".csv")
    maxent_file=read.csv(file_name)
    maxent_file$metric=eval_stat
    assign(paste0("maxent_", eval_stat), maxent_file)
    file_name=paste0("combined_results/model_eval_metric/eval_metric_comparison_GBM_", eval_stat, ".csv")
    gbm_file=read.csv(file_name)
    gbm_file$metric=eval_stat
    assign(paste0("gbm_", eval_stat), gbm_file)
  }
  all_maxent_metrics=rbind(maxent_TSS, maxent_ROC, maxent_KAPPA)
  all_gbm_metrics=rbind(gbm_TSS, gbm_ROC, gbm_KAPPA)
  
  all_maxent_metrics$Species=gsub(pattern="_", replacement=" ", all_maxent_metrics$Species)
  all_gbm_metrics$Species=gsub(pattern="_", replacement=" ", all_gbm_metrics$Species)
  
  a=ggplot(data=all_maxent_metrics, aes(x=Species, y=Regional, fill=metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.title=element_blank()) + ylab("Regional maxent eval. metrics") +xlab("")
  a
  tiff_name=paste0("combined_results/model_eval_metric/regional_HI_maxent_eval_metric_comparison.tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  a=ggplot(data=all_maxent_metrics, aes(x=Species, y=Global, fill=metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.title=element_blank()) + ylab("Global maxent eval. metrics") +xlab("")
  a
  tiff_name=paste0("combined_results/model_eval_metric/global_notHI_maxent_eval_metric_comparison.tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  a=ggplot(data=all_gbm_metrics, aes(x=Species, y=Regional, fill=metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.title=element_blank()) + ylab("Regional gbm eval. metrics") +xlab("")
  a
  tiff_name=paste0("combined_results/model_eval_metric/regional_HI_gbm_eval_metric_comparison.tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  a=ggplot(data=all_gbm_metrics, aes(x=Species, y=Global, fill=metric)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.title=element_blank()) + ylab("Global gbm eval. metrics") +xlab("")
  a
  tiff_name=paste0("combined_results/model_eval_metric/global_notHI_gbm_eval_metric_comparison.tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
}
