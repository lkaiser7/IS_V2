### intial script collecting species data ###
### download online gbif data collections ###
### invasive species global location data ###

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
# create directory for gbif data
dir.create(gbifDir, showWarnings = FALSE

# load necessary installed package
library("rgbif")
library("plyr")
library("ggplot2")
library("stringr")

# invasive species scientific names (16 only - exclude palmgrass [see below])
sp_nm = c('Clidemia hirta', 'Falcataria moluccana', 'Hedychium gardnerianum', 
          'Lantana camara', 'Leucaena leucocephala', 'Melinis minutiflora', 
          'Miconia calvescens', 'Morella faya', 'Panicum maximum', 'Passiflora tarminiana',
          'Pennisetum clandestinum', 'Pennisetum setaceum', 'Psidium cattleianum',
          'Schinus terebinthifolius', 'Cyathea cooperi', 'Ulex europaeus')
# NOTE: Cyathea cooperi is the species synonym for Sphaeropteris cooperi

for (i in 1:length(sp_nm)) {  # set i = 1 for debugging 
  # determine GBIF key for species based on scientific name
  key<-name_backbone(name = sp_nm[i])$speciesKey
  # search records of data for species based on GBIF key
  sp_count<-occ_count(taxonKey = key, georeferenced = TRUE)
  # combine species count information
  sp_info<-data.frame(sp_nm[i], key, sp_count)
  if (i == 1) {
    all_info<-sp_info
  } else {
    all_info<-rbind(all_info, sp_info)
  }
} # END species information for loop

# collect palmgrass separately due to multiple database records
key<-name_backbone(name = 'Setaria palmifolia', kingdom = 'Plantae', 
                    phylum = 'Magnoliophyta', class = 'Liliopsida', 
                    order = 'Poales', family = 'Poaceae',
                    genus = 'Setaria')$speciesKey
sp_count<-occ_count(taxonKey = key, georeferenced = TRUE)
sp_info<-data.frame(sp_nm.i. = 'Setaria palmifolia', key, sp_count)

# combine palmgrass info with all other info for complete database
all_info<-rbind(all_info, sp_info)

# save final output table of all invasive species information for reference
write.csv(all_info, file = paste0(dataDir, "sp_gbif_count.csv"))
# all_info<-read.csv(paste0(dataDir, "sp_gbif_count.csv"), header = TRUE)

# use all_info file to collect data for all 17 species based on GBIF key
for (n in 1:length(all_info$key)) {  # set n = 1 for debugging
  # set n = 4 for Miconia example data set
  # return data fields for each species from GBIF database for first 500 records
  sp_data<-occ_search(taxonKey = all_info$key[n], 
                      hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                      fields = 'minimal', return = "data")
  
  # collect all species data - only records of 500 returned per inquiry
  # create counter for 'start' input and if statement condition for intervals of 500
  rec<-500
  # loop through database search and collect all species location records
  repeat {
    # collect next group of 500 record points per individual species 
    sp_search2<-occ_search(taxonKey = all_info$key[n], 
                           hasCoordinate = TRUE, hasGeospatialIssue = FALSE,
                           fields = 'minimal', return = "data", start = rec)
    # merge additional searches with inital species data object for complete record
    sp_data<-rbind(sp_data, sp_search2)
    # increase counter interval by 500 for next group of species location data
    rec<-rec + 500
    # repeat until all records of species location data have been collected
    if (all_info$sp_count[n] < rec) {
      break
    } # END if condition to exit repeat loop
  } # END repeat loop of species data search and collection
    
  # isolate species name
  sp_name<-sp_data[1, 1]
  # edit species name file naming
  sp_name<-str_replace_all(sp_name, " ", "_")
  
  # remove last row of species data that is NA (not discrete values)
  last_row<-length(sp_data$key)
  # check to see if last row is not valid data entry
  data_chk<-sp_data[last_row, 1] == sp_data[last_row - 1, 1]
  # if invalid, then remove last column of data frame
  if (data_chk == FALSE) {
    sp_data<-sp_data[1:last_row - 1, ]
  }
  
  # save species GBIF location data
  write.csv(sp_data, file = paste0(gbifDir, sp_name, ".csv"))
  # sign-posting of completed data collection per species
  cat('\n All' last_row 'data points collected and saved for', sp_name, '\n')
} # END species data collection
  
# create loop to map all species location data and save as tiff file
for (m in 1:length(all_info$key)) { # set m = 1 for debugging or individual runs
  # edit species name to match file naming convention
  new_sp<-str_replace_all(all_info$sp_nm.i.[m], " ", "_")
  # open species csv file previously saved above
  new_dat<-read.csv(paste0(gbifDir, new_sp, ".csv"), header = TRUE)
  
  # create plot and save as tiff file
  tiff(filename = paste0(gbifDir, new_sp, "_map.png"),
       res = 300, units = "in", pointsize = 12,
       width = 10, height = 10, compression = "lzw")
  # build gbif map
  new_map<-gbifmap(new_dat)
  # print map
  print(new_map + ggtitle(new_sp))
  # save tiff file 
  dev.off()  
} # END species mapping for loop

############################
### END gbif_data SCRIPT ###
############################