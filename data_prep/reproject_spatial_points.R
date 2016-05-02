### reprojection of spatial data points ###
### changing UTM to Lat/Lon WGS84 datum ###
### for local species location datasets ###

# set root path to source files
rootDir<-"Y:/PICCC_analysis/IS_analysis/"

# set working directory to main analysis folder
setwd(rootDir)

# location of scripts and code
codeDir<-"C:/Users/Lauren/Dropbox/GitHub/IS_V2/"

# location of all data
dataDir<-paste0(rootDir, "data/")

# load all raw species data
all_data<-read.csv(paste0(dataDir, "all_sp_local_occ.csv"), header = T, stringsAsFactors = F)

# combine eastings and northings data
coords<-cbind(all_data$X, all_data$Y)

# create data frame of spatial points
spatial_df<-SpatialPointsDataFrame(coords, all_data[, c("X", "Y", "Species")])

# define class of corrdinate reference system to use
proj4string(spatial_df)<-CRS("+init=epsg:26904") # NAD83 zone 4N (& 5N?)
# from http://spatialreference.org/ref/epsg/

# transform datum and map projection to WGS84
sp_df<-spTransform(spatial_df, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# store data output
all_data<-as.data.frame(sp_df)

# rename column headers
names(all_data)<-c("E", "N", "Species", "decimalLongitude", "decimalLatitude")

# save output data
write.csv(all_data, paste0(dataDir, "sp_loc_occ_wgs84.csv"))

############################
### END reproject SCRIPT ###
############################