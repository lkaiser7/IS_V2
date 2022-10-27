### script to format species data at a regional scale ###
### parse data for every species from a given file ###
### all species location data for Hawaiian Islands ###

#load necessary pacakges
library("stringr")

# set root path to source files
rootDir<-"Y:/PICCC_analysis/IS_analysis/"
# set working directory to main analysis folder
setwd(rootDir)

# location of scripts and code
codeDir<-"C:/Users/Lauren/Dropbox/GitHub/IS_V2/"

# location of all data 
dataDir<-paste0(rootDir, "data/")
# regional data
outDir<-paste0(dataDir, "raw_data/regional_data/")
# create directory for regional data
dir.create(outDir, showWarnings = FALSE)

# open raw data file with all species presences location data
all_pres<-read.csv(paste0(dataDir, "sp_loc_occ_wgs84.csv"), header = TRUE)

# separate and parse raw data per species and save output file
sp_split<-split(all_pres, all_pres$Species)
# loop through list and save each list element separately per station
for (i in 1:length(sp_split)) {  # set i = 1 for debugging 
  sp_name<-names(sp_split)[[i]]
  sp_name<-str_replace_all(sp_name, " ","_")
  write.csv(sp_split[[i]], file = paste0(outDir, sp_name, ".csv"))
}

# NOTE: manually replace Spaeropteris with Cyathea (synonym species) file name
# NOTE: manually select Passiflora mollisima data due to naming convention

# loop through new files and determine number of counts per species
sp_files<-list.files(outDir)
for (n in 1:length(sp_files)) {  # set n = 2 for debugging (skip raw data file) 
  # open csv file for each individual species 
  sp_csv<-read.csv(paste0(outDir, sp_files[n]), header = T)
  # determine number of observations per species 
  sp_rec<-length(sp_csv$Species)
  # create species record information
  sp_nm.i.<-sp_csv$Species[1]
  sp_dat<-data.frame(sp_nm.i., sp_rec)
  
  # combine all species data information 
  if (n == 1){
    all_dat<-sp_dat
  } else {
    all_dat<-rbind(all_dat, sp_dat)
  }
}

# save output file of species records 
write.csv(all_dat, file = paste0(dataDir, "sp_regional_count.csv"))

# format BISON data downloaded from: http://bison.usgs.ornl.gov/#home
bisonDir<-paste0(dataDir, "raw_data/bison_data/")
# loop through bison files to edit and determine number of counts per species
b_files<-list.files(bisonDir)
for (j in 1:length(b_files)) {  # set j = 1 for debugging
  # open csv file for each individual species 
  b_csv<-read.csv(paste0(bisonDir, b_files[j]), header = F, na.strings = c("", "NA"))
  # retain only rows and columns needed
  b_csv<-cbind(b_csv[-1, c(14, 21, 20)])
  # keep rows with complete cases (no NAs)
  b_csv<-b_csv[complete.cases(b_csv), ]
  
  # rename headers
  names(b_csv)<-c("Species", "decimalLongitude", "decimalLatitude")
  # define species name
  sp_nm<-b_csv[1, 1]
  sp_nm<-str_replace_all(sp_nm, " ","_")
  
  # save formatted csv file
  write.csv(b_csv, file = paste0(bisonDir, sp_nm, ".csv"))
  
  # determine number of observations per species 
  b_rec<-length(b_csv$Species)
  # create species record information
  sp_nm.i.<-b_csv$Species[1]
  b_dat<-data.frame(sp_nm.i., b_rec)
  
  # combine all species data information 
  if (j == 1){
    b_all_dat<-b_dat
  } else {
    b_all_dat<-rbind(b_all_dat, b_dat)
  }
}

# merge with GBIF data count records
gbif_dat<-read.csv("data/sp_gbif_count.csv", header = TRUE)
count_dat<-merge(gbif_dat, all_dat, by = "sp_nm.i.")
names(count_dat)<-c("Species", "X", "Key", "GBIF", "Regional")

# save combinde output data count
write.csv(count_dat, "data/sp_count_rec.csv")

#############################
### END regional_data SCRIPT ###
#############################