wd="D:/projects/Invasives_modeling/results/xfirst round of results/main_results/"
setwd(wd)

eval_metrics=c("TSS", "ROC", "KAPPA") #

eval_metric = eval_metrics[1]
for (eval_metric in eval_metrics){
  global_tss_df=read.csv(paste0("global_models/outputs/all_eval_mat_",eval_metric,".csv"))
  local_tss_df=read.csv(paste0("local_models/outputs/all_eval_mat_",eval_metric,".csv"))
  
  global_LRG_tss_df=read.csv(paste0("global_models_LRG/outputs/all_eval_mat_",eval_metric,".csv"))
  local_LRG_tss_df=read.csv(paste0("local_models_LRG/outputs/all_eval_mat_",eval_metric,".csv"))
  
  calculate_eval_mean_eval_metric=function(metric_df, model_eval_metric="MAXENT.Phillips"){
    metric_df_maxent_mean=metric_df[metric_df$rownames.Spp_eval.==model_eval_metric, ]
    metric_df_maxent_mean=metric_df_maxent_mean[,-2]
    jnk=apply(metric_df_maxent_mean[,-1], 1, mean, na.rm=T)
    mean_eval_metric=data.frame(species=metric_df_maxent_mean[, 1], eval_metric=jnk)
    return(mean_eval_metric)
  }
  
  mean_global_maxent_eval=calculate_eval_mean_eval_metric(global_tss_df)
  mean_local_maxent_eval=calculate_eval_mean_eval_metric(local_tss_df)
  mean_global_LRG_maxent_eval=calculate_eval_mean_eval_metric(global_LRG_tss_df)
  mean_local_LRG_maxent_eval=calculate_eval_mean_eval_metric(local_LRG_tss_df)
  
  mean_global_GBM_eval=calculate_eval_mean_eval_metric(global_tss_df, model_eval_metric="GBM")
  mean_local_GBM_eval=calculate_eval_mean_eval_metric(local_tss_df, model_eval_metric="GBM")
  mean_global_LRG_GBM_eval=calculate_eval_mean_eval_metric(global_LRG_tss_df, model_eval_metric="GBM")
  mean_local_LRG_GBM_eval=calculate_eval_mean_eval_metric(local_LRG_tss_df, model_eval_metric="GBM")
  
  mean_global_GBM_eval=rbind(mean_global_GBM_eval, mean_global_LRG_GBM_eval)
  mean_local_GBM_eval=rbind(mean_local_GBM_eval, mean_local_LRG_GBM_eval)
  
  mean_global_maxent_eval=rbind(mean_global_maxent_eval, mean_global_LRG_maxent_eval)
  mean_local_maxent_eval=rbind(mean_local_maxent_eval, mean_local_LRG_maxent_eval)
  
  #View(mean_global_LRG_GBM_eval)
  #View(mean_global_GBM_eval)
  
  mean_maxent_eval_df=merge(mean_global_maxent_eval, mean_local_maxent_eval, by="species")
  mean_GBM_eval_df=merge(mean_global_GBM_eval, mean_local_GBM_eval, by="species")
  
  names(mean_maxent_eval_df)=c("Species", "Global", "Local")
  names(mean_GBM_eval_df)=c("Species", "Global", "Local")
  
  file_name=paste0("model_eval_metric/eval_metric_comparison_maxent_", eval_metric, ".csv")
  write.csv(mean_maxent_eval_df, file_name, row.names = F)
  file_name=paste0("model_eval_metric/eval_metric_comparison_GBM_", eval_metric, ".csv")
  write.csv(mean_GBM_eval_df, file_name, row.names = F)
  
  plot(mean_maxent_eval_df$Global, mean_maxent_eval_df$Local)
  plot(mean_GBM_eval_df$Global, mean_GBM_eval_df$Local)
  
  geom.text.size = 2
  theme.size = (14/5) * geom.text.size
  
  library(ggplot2)
  a=ggplot(mean_maxent_eval_df, aes(x=Global, y=Local)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=mean_maxent_eval_df$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")
  a
  tiff_name=paste0("model_eval_metric/eval_metric_comparison_maxent_", eval_metric, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
  library(ggplot2)
  a=ggplot(mean_GBM_eval_df, aes(x=Global, y=Local)) + 
    geom_point(aes(size=1.25)) +
    geom_text(label=mean_GBM_eval_df$Species, nudge_x = 0.0, nudge_y = 0.015,  size=geom.text.size)+ 
    theme(legend.position="none")
  a
  tiff_name=paste0(wd, "model_eval_metric/eval_metric_comparison_GBM_", eval_metric, ".tiff")
  ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
  
}

#now compare metrics
for (eval_metric in eval_metrics){
  file_name=paste0("model_eval_metric/eval_metric_comparison_maxent_", eval_metric, ".csv")
  maxent_file=read.csv(file_name)
  maxent_file$metric=eval_metric
  assign(paste0("maxent_", eval_metric), maxent_file)
  file_name=paste0("model_eval_metric/eval_metric_comparison_GBM_", eval_metric, ".csv")
  gbm_file=read.csv(file_name)
  gbm_file$metric=eval_metric
  assign(paste0("gbm_", eval_metric), gbm_file)
}
all_maxent_metrics=rbind(maxent_TSS, maxent_ROC, maxent_KAPPA)
all_gbm_metrics=rbind(gbm_TSS, gbm_ROC, gbm_KAPPA)

all_maxent_metrics$Species=gsub(pattern="_", replacement=" ", all_maxent_metrics$Species)
all_gbm_metrics$Species=gsub(pattern="_", replacement=" ", all_gbm_metrics$Species)

a=ggplot(data=all_maxent_metrics, aes(x=Species, y=Local, fill=metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title=element_blank()) + ylab("Local maxent eval. metrics") +xlab("")
a
tiff_name=paste0(wd, "model_eval_metric/local_maxent_eval_metric_comparison.tiff")
ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")

a=ggplot(data=all_maxent_metrics, aes(x=Species, y=Global, fill=metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title=element_blank()) + ylab("Global maxent eval. metrics") +xlab("")
a
tiff_name=paste0(wd, "model_eval_metric/global_maxent_eval_metric_comparison.tiff")
ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")

a=ggplot(data=all_gbm_metrics, aes(x=Species, y=Local, fill=metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title=element_blank()) + ylab("Local gbm eval. metrics") +xlab("")
a
tiff_name=paste0(wd, "model_eval_metric/local_gbm_eval_metric_comparison.tiff")
ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")

a=ggplot(data=all_gbm_metrics, aes(x=Species, y=Global, fill=metric)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title=element_blank()) + ylab("Global gbm eval. metrics") +xlab("")
a
tiff_name=paste0(wd, "model_eval_metric/global_gbm_eval_metric_comparison.tiff")
ggsave(filename = tiff_name, plot = a, width = 6, height = 4, units = "in", compress="lzw")
