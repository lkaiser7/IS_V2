### species records count ###
### variable importance by model ###

# clear the environment, temp files, and all other variables
#rm(list = ls())

##################
##### SET UP #####
##################

# set root path to source files
#rootDir<-"D:/Phase1_SDMs/"
# set working directory to main analysis folder
setwd(rootDir)

# set path to data files
#dataDir<-paste0(rootDir, "data/")

#########################
##### SPECIES COUNT #####
#########################

library(rgdal)
# load location coast shapefile from file path to personal database
loc_coast<-readOGR("data/map_data", "Main_Hawaiian_Islands_simple3")
loc_coast<-loc_coast[which(loc_coast$Island != "NI"),]
plot(loc_coast)
graphics.off()

# list species files
sp_files<-list.files(paste0(dataDir, "merged_data/all_data/"))
#sp_files<-sp_files[-8]

# loop through species and table counts

s=1
for(s in 1:length(sp_files)){  # s = 1 for debugging
  cat("plotting HI spp points for species ", s, "\n")
  # open species file
  sp_rec<-read.csv(paste0(dataDir, "merged_data/all_data/", sp_files[s]))
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

#View(all_sp_count)
# save species count file
write.csv(all_sp_count, paste0(dataDir, "all_sp_count.csv"))
rm(s)

###############################
##### VARIABLE IMPORTANCE #####
###############################

# list model scales
all_mod_scale = c("global_notHI", "local_HI", "nested_HI")#c("local", "global", "nested")
all_mod_scale=paste0(all_mod_scale, "_models")

#all_mod_scale = c("local", "global", "nested")

# list bioclims selected
bc_list = gsub("\\..*", "", env_var_files) #c("bio1", "bio7", "bio12", "bio15")


# loop through each model type 
m=1
for(m in 1:length(all_mod_scale)){ # set m = 1 for debugging 
  cat("summarizing var imp for ", all_mod_scale[m], "\n")
  # set model scale
  mod_scale<-all_mod_scale[m]
  
  # set path of ongoing project run for all outputs
  project_path<-paste0(rootDir, mod_scale, "/")
  
  # location to save any extra or more specific outputs
  outDir<-paste0(project_path, "outputs/")
  
  # open variable importance file
  var_imp<-read.csv(paste0(outDir, "all_VariImp.csv"))
  # head(var_imp)
  #View(var_imp)
  
  # loop through each species
  s=1
  for(s in 1:length(all_sp_nm)){ # set s = 1 for debugging
    cat("doing species ", all_sp_nm[s], "\n")
    
    # select single species name
    sp.name<-all_sp_nm[s]
    
    # select species variable importance data
    sp_var_imp<-var_imp[which(var_imp[,1] == sp.name),]
    
    # # separate maxent and gbm runs
    # sp_var_max<-sp_var_imp[grep("MAXENT", names(sp_var_imp))]
    # sp_var_gbm<-sp_var_imp[grep("GBM", names(sp_var_imp))]
    for (model_type in models_to_run){
      sp_var_mod<-sp_var_imp[grep(model_type, names(sp_var_imp))]
      sp_var_mod= rowMeans(sp_var_mod)
      if (model_type == models_to_run[1]){
        sp_bc_vars<-data.frame(SPECIES = sp.name, VAR = bc_list, 
                               sp_var_mod)
      }else{
        sp_bc_vars=cbind(sp_bc_vars, sp_var_mod)
      }
    }
    
    names(sp_bc_vars)=c("SPECIES", "VAR", paste0(models_to_run, "_MEAN"))
    
    # calculate means per variable
    # sp_bc_vars<-cbind(sp.name, bc_list, rowMeans(sp_var_max), rowMeans(sp_var_gbm))
    # sp_bc_vars<-data.frame(SPECIES = sp.name, VAR = bc_list, 
    #                        MAXENT_MEAN = rowMeans(sp_var_max), 
    #                        GBM_MEAN = rowMeans(sp_var_gbm))
    if(s == 1){
      all_bc_vars<-sp_bc_vars
    }else{
      all_bc_vars<-rbind(all_bc_vars, sp_bc_vars)
    }
    
  } # END s loop
  #View(all_bc_vars)
  # save means by model runs
  write.csv(all_bc_vars, paste0(outDir, "all_VariImp_model_mean.csv"))
  
} # END m loop
rm(s)



### END SPECIES COUNTS AND VARIABLE IMPORTANCE ###