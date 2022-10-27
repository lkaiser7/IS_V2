rm(list = ls()) # clear the environment, temp files, and all other variables
# set root path to source files
project_dirs=c("C:/Users/lkaiser-regional/Desktop/Phase1_SDMs/", "E:/invasives_SDM5/", 
               "~/projects/invasives_SDM5/", "/home/pierc/projects/invasives_SDM/")
rootDir=project_dirs[min(which(dir.exists(project_dirs)))]
setwd(rootDir) # set working directory to main analysis folder

# # location of scripts and code
# codeDirs=c("D:/projects/Invasives_modeling/IS_V2_repo/", paste0(rootDir, "IS_V2/"), "/home/pierc/git_repos/IS_V2/") #in order of priority
# codeDir=codeDirs[min(which(dir.exists(codeDirs)))]

# location of all data
dataDir<-paste0(rootDir, "data/")

# # location of map data and shapefiles
# mapDirs<-c(paste0(dataDir, "map_data/"), "D:/data/")
# mapDir<-mapDirs[min(which(dir.exists(mapDirs)))]

bioclims_dirs=c("D:/data/global_climate/wc2.1_30s_bio_simplified/", "/home/pierc/data/global_climate/wc2.1_30s_bio_simplified/", paste0(dataDir, "bioclim_vars/")) #in order of priority
global_bioclims_dir<-bioclims_dirs[min(which(dir.exists(bioclims_dirs)))]
bioclims_dirs=c("/home/pierc/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/", 
                  "D:/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/", 
                  paste0(dataDir, "bioclim_vars/")) #in order of priority
regional_bioclims_dir<-bioclims_dirs[min(which(dir.exists(bioclims_dirs)))]

fitting_bios_global<-paste0(global_bioclims_dir, "all_baseline/current_30s/")
fitting_bios_HIs<-c(paste0(regional_bioclims_dir, "all_HRCM/current_250m_redone/"), "D:/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/",
                    "/home/pierc/data/climate_data/20201123_HRCM_NCAR_projections2/bioclims/baseline_rasters/")
fitting_bios_HI<-fitting_bios_HIs[min(which(dir.exists(fitting_bios_HIs)))]


# select current data and bioclims to use for model approach
global_biofitRun<-fitting_bios_global     # for model fitting
global_bioclim_scaling_factors=c(1/10, 1/10, 1, 1/10) #rasters were saved as integers
regional_biofitRun<-fitting_bios_HI     # for model fitting
regional_bioclim_scaling_factors=c(1, 1, 1, 1)

library("terra")

env_var_files = c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
vars_to_include=c(T, F, T, F)
#vars_to_include=c(T, T, T, T)

env_var_file = env_var_files[1]
j=1
for (env_var_file in env_var_files){
  var_to_include=vars_to_include[j]
  if (var_to_include){
    cat("doing ", env_var_file, "\n")
    global_predictor = rast(paste0(global_biofitRun, env_var_file)) 
    global_predictor=global_predictor*global_bioclim_scaling_factors[j]
    
    regional_predictor = rast(paste0(regional_biofitRun, env_var_file)) 
    regional_predictor=regional_predictor*regional_bioclim_scaling_factors[j]
    g_min=c(global(global_predictor, min, na.rm=T))[[1]]
    g_max=global(global_predictor, max, na.rm=T)[[1]]
    l_min=global(regional_predictor, min, na.rm=T)[[1]]
    l_max=global(regional_predictor, max, na.rm=T)[[1]]
    if (l_min>g_min) {
      min_seq=seq(l_min, g_min, length.out=10)
    }else{
      min_seq=seq(l_min, l_min, length.out=10)
    }
    if (l_max<g_max) {
      max_seq=seq(l_max, g_max, length.out=10)
    }else{
      max_seq=seq(l_max, l_max, length.out=10)
    }
    
    total_seq=c(rev(min_seq), max_seq) 
    total_seq[1]=total_seq[1]-0.1
    total_seq[length(total_seq)]=total_seq[length(total_seq)]+0.1
    class_DF=data.frame(from=total_seq[c(1:(length(total_seq)-1))], to=total_seq[c(2:length(total_seq))], becomes=c(1:10,(10-1):1)/10)
    
    global_predictor_sre=classify(global_predictor, as.matrix(class_DF))
    if (j==1){
      global_predictor_sre_stack=global_predictor_sre
    }else{
      global_predictor_sre_stack = c(global_predictor_sre_stack, global_predictor_sre)
    }    
    
  }
  j=j+1
}

SRE_min=min(global_predictor_sre_stack)
plot(SRE_min)
plot(SRE_min==1)
writeRaster(SRE_min*10, paste0(global_bioclims_dir, "Hawaii_SRE_bio1_12.tif"), overwrite=TRUE, gdal=c("COMPRESS=LZW"), datatype='INT1U')


