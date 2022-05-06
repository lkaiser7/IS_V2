### species records count ###
### variable importance by model ###

# clear the environment, temp files, and all other variables
rm(list = ls())

##################
##### SET UP #####
##################

# set root path to source files
rootDir<-"D:/Phase1_SDMs/"
# set working directory to main analysis folder
setwd(rootDir)

# set path to data files
dataDir<-paste0(rootDir, "data/")

#########################
##### SPECIES COUNT #####
#########################

library(rgdal)
# load location coast shapefile from file path to personal database
loc_coast<-readOGR("data/map_data", "Main_Hawaiian_Islands_simple3")
loc_coast<-loc_coast[which(loc_coast$Island != "NI"),]
plot(loc_coast)

# list species files
sp_files<-list.files(paste0(dataDir, "all_data/"))
sp_files<-sp_files[-8]

# loop through species and table counts
for(s in 1:length(sp_files)){  # s = 1 for debugging
  # open species file
  sp_rec<-read.csv(paste0(dataDir, "all_data/", sp_files[s]))
  # head(sp_rec)

  # table species data
  sp_count<-c(sp_rec$Species[1], table(sp_rec$dataset), dim(sp_rec)[1])

  # create final data set
  if(s == 1){
    all_sp_count<-sp_count
    plot(loc_coast)
    points(sp_rec$decimalLongitude, sp_rec$decimalLatitude, pch = 20, col = s)
  }else{
    all_sp_count<-rbind(all_sp_count, sp_count)
    points(sp_rec$decimalLongitude, sp_rec$decimalLatitude, pch = 20, col = s)
  }

}

# save species count file
write.csv(all_sp_count, paste0(dataDir, "all_sp_count.csv"))
rm(s)

###############################
##### VARIABLE IMPORTANCE #####
###############################

# list model scales
all_mod_scale = c("local", "global", "nested")

# list bioclims selected
bc_list = c("bio1", "bio7", "bio12", "bio15")

# set all_sp_nm = 'Clidemia_hirta' for testing and debugging
# list all species names to be analyzed
sm_sp_nm = c('Clidemia_hirta', 'Falcataria_moluccana', 'Hedychium_gardnerianum',
             'Lantana_camara', 'Leucaena_leucocephala', 'Melinis_minutiflora',
             'Morella_faya', 'Panicum_maximum',
             'Passiflora_tarminiana', 'Pennisetum_clandestinum', 'Pennisetum_setaceum',
             'Psidium_cattleianum', 'Setaria_palmifolia','Schinus_terebinthifolius')

lrg_sp_nm = c('Cyathea_cooperi', 'Miconia_calvescens', 'Ulex_europaeus')
lrg_sp_nm = c('Cyathea_cooperi', 'Ulex_europaeus')

# select current data and bioclims to use for model approach
all_sp_nm<-sm_sp_nm      # small or large species data files 

# loop through each model type 
for(m in 1:length(all_mod_scale)){ # set m = 1 for debugging 
  # set model scale
  mod_scale<-all_mod_scale[m]
  
  # select name for project and create directory
  if(mod_scale == "global"){
    project_run<-paste0(mod_scale, "_notHI_models")
    #project_run<-paste0(mod_scale, "_notHI_models_LRG")
  }else{
    project_run<-paste0(mod_scale, "_HI_models")
    #project_run<-paste0(mod_scale, "_HI_models_LRG")
  }
  # "global_notHI_models", "global_notHI_models_LRG"
  # "local_HI_models", "local_HI_models_LRG"
  # "nested_HI_models", "nested_HI_models_LRG"
  
  # set path of ongoing project run for all outputs
  project_path<-paste0(rootDir, project_run, "/")
  
  # location to save any extra or more specific outputs
  outDir<-paste0(project_path, "outputs/")
  
  # open variable importance file
  var_imp<-read.csv(paste0(outDir, "all_VariImp.csv"))
  # head(var_imp)
  
  # loop through each species
  for(s in 1:length(all_sp_nm)){ # set s = 1 for debugging
    # select single species name
    sp.name<-all_sp_nm[s]
    
    # select species variable importance data
    sp_var_imp<-var_imp[which(var_imp[,1] == sp.name),]
    
    # separate maxent and gbm runs
    sp_var_max<-sp_var_imp[grep("MAXENT", names(sp_var_imp))]
    sp_var_gbm<-sp_var_imp[grep("GBM", names(sp_var_imp))]
    
    # calculate means per variable
    # sp_bc_vars<-cbind(sp.name, bc_list, rowMeans(sp_var_max), rowMeans(sp_var_gbm))
    sp_bc_vars<-data.frame(SPECIES = sp.name, VAR = bc_list, 
                           MAXENT_MEAN = rowMeans(sp_var_max), 
                           GBM_MEAN = rowMeans(sp_var_gbm))
    if(s == 1){
      all_bc_vars<-sp_bc_vars
    }else{
      all_bc_vars<-rbind(all_bc_vars, sp_bc_vars)
    }
    
  } # END s loop
  
  # save means by model runs
  write.csv(all_bc_vars, paste0(outDir, "all_VariImp_model_mean.csv"))
  
} # END m loop
rm(s)

# load necessary packages
library(ggplot2)
library(cowplot)

# create folder in output directory for plots
vi_fold<-paste0(rootDir, "mean_VariImp_plots/")
dir.create(vi_fold, showWarnings = FALSE)

# LOAD FILES
global_vi<-read.csv("D:/Phase1_SDMs/global_notHI_models/outputs/all_VariImp_model_mean.csv")
#global_vi<-read.csv("D:/Phase1_SDMs/global_notHI_models_LRG/outputs/all_VariImp_model_mean.csv")
local_vi<-read.csv("D:/Phase1_SDMs/local_HI_models/outputs/all_VariImp_model_mean.csv")
#local_vi<-read.csv("D:/Phase1_SDMs/local_HI_models_LRG/outputs/all_VariImp_model_mean.csv")
nested_vi<-read.csv("D:/Phase1_SDMs/nested_HI_models/outputs/all_VariImp_model_mean.csv")
#nested_vi<-read.csv("D:/Phase1_SDMs/nested_HI_models_LRG/outputs/all_VariImp_model_mean.csv")

# select current data and bioclims to use for model approach
all_sp_nm<-sm_sp_nm      # small or large species data files 
#all_sp_nm<-lrg_sp_nm      # small or large species data files 

# loop through each species and format data to plot
for(s in 1:length(all_sp_nm)){
  # select species
  sp_nm<-all_sp_nm[s]
  
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

### END SPECIES COUNTS AND VARIABLE IMPORTANCE ###