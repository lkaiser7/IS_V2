### combining all collected species data ###
### gbif, local & bison data collections ###
### species location data to be analyzed ###

# set root path to source files
rootDir<-"Y:/PICCC_analysis/IS_analysis/"
# set working directory to main analysis folder
setwd(rootDir)

# location of scripts and code
codeDir<-"C:/Users/Lauren/Dropbox/GitHub/IS_V2/"

# location of all data
dataDir<-paste0(rootDir, "data/")

# global GBIF data
gbifDir<-paste0(dataDir, "raw_data/gbif_data/")
# regional BISON data
bisonDir<-paste0(dataDir, "raw_data/bison_data/")
# local internal data
localDir<-paste0(dataDir, "raw_data/local_data/")

# output folder for all combinded data 
allDir<-paste0(dataDir, "all_data/")
# output folder for Hawaii data
hiDir<-paste0(dataDir, "hi_data/")
# output folder excluding Hawaii data
nohiDir<-paste0(dataDir, "no_hi_data/")

# list all invasive species 
all_sp_nm = c('Clidemia_hirta', 'Falcataria_moluccana', 'Hedychium_gardnerianum', 
              'Lantana_camara', 'Leucaena_leucocephala', 'Melinis_minutiflora', 
              'Miconia_calvescens', 'Morella_faya', 'Panicum_maximum', 
              'Passiflora_tarminiana', 'Pennisetum_clandestinum', 'Pennisetum_setaceum', 
              'Psidium_cattleianum', 'Setaria_palmifolia','Schinus_terebinthifolius', 
              'Cyathea_cooperi', 'Ulex_europaeus')
# NOTE: Cyathea cooperi is the species synonym for Sphaeropteris cooperi
# NOTE: Passiflora tarminiana is a species synonym of Passiflora mollisima

# loop through each file per species 
for(k in 1:length(all_sp_nm)) {  # set k = 1 for debugging
  # open global species data file
  g_data<-read.csv(paste0(gbifDir, all_sp_nm[k], ".csv"), header = TRUE)
  names(g_data)[1]<-"Species"  # column name formatting
  # open bison species data file
  b_data<-read.csv(paste0(bisonDir, all_sp_nm[k], ".csv"), header = TRUE)
  # open local species data file
  l_data<-read.csv(paste0(localDir, all_sp_nm[k], ".csv"), header = TRUE)
  
  # keep selected columns of data from each dataset
  g_data<-with(g_data, { # Species[2], decimalLatitude[4], decimalLongitude[5]
    data.frame(Species, decimalLongitude, decimalLatitude, rep("G", dim(g_data)[1]))})
  names(g_data)[4]<-"dataset"  # column name formatting
  b_data<-with(b_data, { # Species[2], decimalLongitude[3], decimalLatitude[4]
    data.frame(Species, decimalLongitude, decimalLatitude, rep("B", dim(b_data)[1]))})
  names(b_data)[4]<-"dataset"  # column name formatting
  l_data<-with(l_data, { # Species[5], decimalLongitude[6], decimalLatitude[7]
    data.frame(Species, decimalLongitude, decimalLatitude, rep("L", dim(l_data)[1]))})
  names(l_data)[4]<-"dataset"  # column name formatting
  
  # merge data into single data frame per species
  all_sp_data<-rbind(g_data, b_data, l_data)
  # replace species column with name for consistency
  all_sp_data$Species<-rep(all_sp_nm[k], dim(all_sp_data)[1])
  
  # save combined species data
  write.csv(all_sp_data, file = paste0(allDir, all_sp_nm[k], ".csv"))
  
  # limit the longitude (-160 < X < -150) & latitude (15 < Y < 25) extent
  hi_sp_data<-subset(all_sp_data, decimalLongitude > -160 & decimalLongitude < -150 
                     & decimalLatitude > 15 & decimalLatitude < 25)
  # save Hawaii data
  write.csv(hi_sp_data, file = paste0(hiDir, all_sp_nm[k], ".csv"))
  
  # exlude all Hawaii data 
  nohi_data<-subset(all_sp_data, decimalLongitude < -160 | decimalLongitude > -150 
                     | decimalLatitude < 15 | decimalLatitude > 25)
  # save data exluding Hawaii
  write.csv(nohi_data, file = paste0(nohiDir, all_sp_nm[k], ".csv"))
  
} 
# end loop through species

#############################
### END merge_data SCRIPT ###
#############################