rm(list = ls()) # clear the environment, temp files, and all other variables
seed=9214621 #NULL
cpucores = 9 # select number of computer cores for processing (max = 32)
eval_stats = c("TSS") #DEBUG
run_type="regional_HI" # global_notHI regional_HI nested_HI # select name for project and create directory
run_scripts_beyond_projection=F
run_scripts_after_3_model_scales_done=F

all_sp_nm = c('Clidemia_hirta', 'Falcataria_moluccana', 'Hedychium_gardnerianum',
              'Lantana_camara', 'Leucaena_leucocephala', 'Melinis_minutiflora',
              'Morella_faya', 'Panicum_maximum',
              'Passiflora_tarminiana', 'Pennisetum_clandestinum', 'Pennisetum_setaceum',
              'Psidium_cattleianum', 'Setaria_palmifolia','Schinus_terebinthifolius',
              'Cyathea_cooperi', 'Miconia_calvescens', 'Ulex_europaeus')
project_dirs=c("C:/Users/lkaiser-regional/Desktop/Phase1_SDMs/", "E:/invasives_SDM10/",  "D:/projects/invasives_SDM10/", 
               "/hdd/projects/invasives_SDM10/")
rootDir=project_dirs[min(which(dir.exists(project_dirs)))]
setwd(rootDir) # set working directory to main analysis folder
project_run<-paste0(run_type, "_models")
# set path of ongoing project run for all outputs
project_path<-paste0(rootDir, project_run, "/")
library(stringr)

# sp_nm=all_sp_nm[1]
for (sp_nm in all_sp_nm){
  sp_nm = as.character(sp_nm) 
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  
  SP_ALL_data=read.csv(paste0(project_path, sp_dir, sp_nm, "_bioclim_points.csv"))
  SP_ALL_data=SP_ALL_data[SP_ALL_data$PA %in% c(1),]
  SP_ALL_data=cbind(sp_nm, SP_ALL_data)
  
  if (sp_nm==all_sp_nm[1]){
    SPP_ALL_data=SP_ALL_data
  }else{
    SPP_ALL_data=rbind(SPP_ALL_data, SP_ALL_data)
  }
}

#View(SPP_ALL_data)
library(terra)
DEM_30m=rast("D:/data/DEM_data_HI/30m_2023/DEM_30m_rounded_WGS84_geotif.tif")
# DEM_30m_wgs84 <- project(DEM_30m, "+proj=longlat +datum=WGS84", method="near")
# writeRaster(DEM_30m_wgs84, "D:/data/DEM_data_HI/30m_2023/DEM_30m_rounded_WGS84_geotif.tif", overwrite=TRUE, gdal=c("compress=LZW"), datatype="INT2U")
tmp<-terra::extract(DEM_30m, SPP_ALL_data[, c("X","Y")], cells = F) 
SPP_ALL_data$elev=tmp$Layer_1

library(ggplot2)
df=SPP_ALL_data
ggplot(df, aes(x = elev)) + 
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  facet_wrap(~ sp_nm) + 
  theme_minimal() +
  labs(x = "Elevation", y = "Frequency", title = "Elevation Histograms by Species") +
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))


