##############################
##############################
#run this only after running 1c for all global, local and nested models

# load necessary packages
library(ggplot2)
library(cowplot)

# create folder in output directory for plots
vi_fold<-paste0(rootDir, "combined_results/mean_VariImp_plots/")
dir.create(vi_fold, showWarnings = FALSE, recursive = T)

# LOAD FILES
global_vi<-read.csv("global_notHI_models/outputs/all_VariImp_model_mean.csv")
local_vi<-read.csv("local_HI_models/outputs/all_VariImp_model_mean.csv")
nested_vi<-read.csv("nested_HI_models/outputs/all_VariImp_model_mean.csv")

# loop through each species and format data to plot
s=1
for(s in 1:length(all_sp_nm)){
  # select species
  sp_nm<-all_sp_nm[s]
  cat("doing ", sp_nm, "\n")
  
  # select species data
  sp_data<-rbind(global_vi[which(global_vi$SPECIES == sp_nm),],
                 local_vi[which(local_vi$SPECIES == sp_nm),],
                 nested_vi[which(nested_vi$SPECIES == sp_nm),])
  # add column for model runs
  sp_data$MODEL_RUN<-c(rep("Global", 4), rep("Local", 4), rep("Nested", 4))
  
  # color blind palette
  cb_palette<-c("#009E73", "#E69F00", "#CC79a7", "#56B4E9")
  
  # plot variable importance by model run per species
  # MAXENT Plot
  maxent_vi<-ggplot(sp_data, aes(VAR, MAXENT_MEAN, fill = VAR)) + 
    geom_col() + facet_wrap(~ MODEL_RUN) + ylim(0, 1) +
    scale_fill_manual(values = cb_palette) + 
    labs(x = NULL, y = 'Mean Variable Importance') + ggtitle(paste(sp_nm, "MaxEnt Models")) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  # GBM Plot
  gbm_vi<-ggplot(sp_data, aes(VAR, GBM_MEAN, fill = VAR)) + 
    geom_col() + facet_wrap(~ MODEL_RUN) + ylim(0, 1) +
    scale_fill_manual(values = cb_palette) + 
    labs(x = NULL, y = 'Mean Variable Importance') + ggtitle(paste(sp_nm, "GBM Models")) + 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  
  # arrange plots and save
  plot_grid(maxent_vi, gbm_vi, ncol = 1)
  ggsave(filename = paste0(vi_fold, sp_nm, "_mean_variable_importance_bar_plot.tiff"))
}
