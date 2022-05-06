### BIOMOD2: climate envelop ###
### see source script inputs ###
### step 5:

##################
##### SET UP #####
##################

# set working directory
setwd(project_path)
# load necessary packages 
library(biomod2)
library(stringr)

# bioclims to use for baseline scenario
if (toCompareWithCurrentClimate == 1){
  clim_surface_to_use = biobaseRun 
  proj_nm0 = 'baseline'}
# bioclims to use for future scenario
if (toCompareWithCurrentClimate == 4){
  clim_surface_to_use = biofutureRun 
  proj_nm0 = 'future'}

# create empty vector for variable names
var_name = c()
# loop through all variable file names
for (env_var_file  in env_var_files){
  # remove file extension
  a = strsplit(env_var_file, "\\.")
  # save variable names in vector
  var_name = c(var_name, a[[1]][1])
}

# set memory.limit
(size = 240000)

# create directory for climate analogs in outputs folder
dir.create("output_rasters/", showWarnings = F)
dir.create("output_rasters/analog_clim/", showWarnings = F)

#####################################
##### CREATE ANALOG CLIMATE MAP #####
#####################################

# create empty vector for number of absences to be removed
n_abs_removed = c()

# print sign posting of ongoing climate envelope mapping
cat('\n future analog climate mapping...')

# store file output name for raster
raster_name = paste0("output_rasters/analog_clim/", 
                     comp_projects[2], "_analog_climate")
# store file output name for tiff image
tiff_name = paste0(raster_name, ".tif")
# store file output name for jpeg image
jpeg_name = paste0(raster_name, ".jpg")

# check to see if the analysis for this species was already done
if (file.exists(tiff_name)==F | overwrite==1){ 
  
  # load vector of climate data names
  all_clim_data = c("biobaseRun", "clim_surface_to_use")
  # load vector of bioclimate variable stack
  clim_stacks = c("biobaseRun", "biofutureRun")
  
  # load bioclimate fitting variables 
  clim_data_dir = biofitRun 
  # create temporary vector length of number of bioclims selected
  jnk0 = length(env_var_files)
  
  # load all baseline and comparative predictor variables
  for (dfd in 1:length(all_clim_data)){
    # load and get climate data for projection run
    clim_data = all_clim_data[dfd]
    clim_data_dir = get(clim_data)
    # load predictors and crop to map extent
    predictors_temp = raster(paste0(clim_data_dir, env_var_files[1]))
    #predictors_temp = crop(predictors_temp,  crop_ext)
    # add additional predictors and crop to same extent
    for (jj in 2:jnk0){
      temp = raster(paste0(clim_data_dir, env_var_files[jj]))
      #temp = crop(temp,  crop_ext)
      # stack predictors togetger in layers
      predictors_temp = addLayer(predictors_temp, temp)
    }
    # rename raster stack with variable names
    names(predictors_temp)<-var_name
    # store predictors for projection run 
    assign(clim_stacks[dfd], predictors_temp)
  }
  
  # store single bioclim variable as raster and crop to extent
  jnk = raster(paste0(clim_data_dir, "bio12.tif"))
  #jnk = crop(jnk,  crop_ext)
  # reclassify raster values of response variable 
  Response_var = reclassify(jnk,  c(0.1,Inf,1))
  
  # run surface range envelop analysis on bioclim variables 
  pred<-sre(Response_var, biobaseRun, biofutureRun, Quant = 0.0001)
  # select comparative variable climate analogs
  analog_climates2100=subset(pred, 1)
  # remove temporary vectors and raster stacks
  rm("predictors_temp", "temp", "jnk")
  
  # create blank jpeg file
  jpeg(jpeg_name, width = 10, height = 10, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  # plot future analog climates
  plot(analog_climates2100, legend = FALSE, 
       main = "Climate Envelop",
       xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  plot(map_to_use, add = TRUE, border = "slategray")
  legend("bottomleft", legend = "Analog Climates",
         pch = 15, col = "darkgreen")    
  # save image file
  dev.off()
  
  # create blank tiff file
  tiff(tiff_name, res = 300, units = "in", pointsize = 12,
       width = 10, height = 10, compression = "lzw")
  # plot future analog climates
  plot(analog_climates2100, legend = FALSE, 
       main = "Climate Envelop",
       xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
  plot(map_to_use, add = TRUE, border = "slategray")
  legend("bottomleft", legend = "Analog Climates",
         pch = 15, col = "darkgreen")  
  # save image file
  dev.off() 
  
  # save raster
  writeRaster(analog_climates2100, paste0(raster_name, "_raster.tif"), 
              format = "GTiff", overwrite = TRUE)
  
  # print sign posting of completed climate envelope mapping
  cat('\n future analog climate raster created.')
}  

######################################
##### END ANALOG CLIMATE MAPPING #####
######################################