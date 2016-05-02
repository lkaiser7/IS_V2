# load packages
library(rgdal)
library(raster)
library(stringr)
library(tools)

# set project path and working directory
proj_path<-"Y:/PICCC_analysis/IS_analysis/data/bioclim_vars/"
setwd(proj_path)

# select scenario path 
base<-paste0(proj_path, "all_baseline/")

# download bioclims from WorldClim (http://www.worldclim.org/download) to path
raw_bios2<-getData('worldclim', var = 'bio', res = 2.5, path = base)
raw_bios5<-getData('worldclim', var = 'bio', res = 5, path = base)
raw_bios10<-getData('worldclim', var = 'bio', res = 10, path = base)
# see ?getData for other variables and inputs (res = 2.5, 5. or 10 minutes of a degree)

# set path to data to be processed 
base30s<-paste0(base, "current_30s/")
base2min<-paste0(base, "current_2min/")
base5min<-paste0(base, "current_5min/")
base10min<-paste0(base, "current_10min/")

# function to convert .bil files to .tif files for modeling
bil_to_tif<-function(current_path) { # set current_path = base30s for debugging
  # load all bioclim variables from specific resolution
  bc_vars<-grep(".bil", list.files(paste0(current_path, "raw/")), value = TRUE)
  # loop through 19 bioclims *WARNING* R READS NUMBER VALUES DIFFERNTLY
  # EX: ORDER WILL BE 1="bio_1.bil", 2="bio_10.bil", 3="bio_11.bil", etc.
  for (bio in 1:length(bc_vars)) { #set bio = 2 for debugging
    # split file name for proper nomaclature to avoid misnaming files 
    bil_name<-bc_vars[bio]
    # save bioclim name for .tif file
    tif_name<-str_replace_all(file_path_sans_ext(bil_name), "_", "")
    # rasterize bioclimatic variables 
    bio_var<-raster(paste0(current_path, "raw/", bc_vars[bio]))
    # save raster as .tif file
    writeRaster(bio_var, filename = paste0(current_path, tif_name, ".tif"), 
                format = "GTiff", overwrite = TRUE)
  }
}
# run function and repeat for each spatial resolution of bioclims
bil_to_tif(base30s)
bil_to_tif(base2min)
bil_to_tif(base5min)
bil_to_tif(base10min)
# lower resolutions take longer processing time to write .tif files

### END PROCESSING OF BASELINE BIOCLIMATIC VARIABLES ### 

# plot raster stacks to check outputs of rasters to make sure they are named properly
check_30s<-stack(paste0(base30s, grep('.tif', list.files(base30s), value = T)))
check_2min<-stack(paste0(base2min, grep('.tif', list.files(base2min), value = T)))
check_5min<-stack(paste0(base5min, grep('.tif', list.files(base5min), value = T)))
check_10min<-stack(paste0(base10min, grep('.tif', list.files(base10min), value = T)))

# set Hawaiian extent
hi_ext<-c(-161, -154, 19, 23)
# crop and plot to Hawaiian extent to check patterns of variables (see HRCM)
hi_30s<-crop(check_30s, hi_ext)
hi_2min<-crop(check_2min, hi_ext)
hi_5min<-crop(check_5min, hi_ext)
hi_10min<-crop(check_10min, hi_ext)

### EDIT HRCM BIOCLIMS TO MAINTAIN CONSISTENT UNITS ###
# path for HRCM bioclims
hrcm<-paste0(proj_path, "all_HRCM/")
# set locations for each scenario
c125<-paste0(hrcm, "current_125m/")
c500<-paste0(hrcm, "current_500m/")
f500<-paste0(hrcm, "future_500m/")

# function to convert decimal units of degrees celsius to integers for consistency 
edit_units_tif<-function(hrcm_path) { # set hrcm_path = c125 for debugging
  # load all bioclim variables from specific resolution
  hrcm_vars<-grep(".tif", list.files(paste0(hrcm_path, "raw/")), value = TRUE)
  # loop through selected bioclims in raw folder with original data 
  for (hrcm_bio in 1:length(hrcm_vars)) { #set hrcm_bio = 1 for debugging
    # load original .tif file as raster
    og_tif<-raster(paste0(hrcm_path, "raw/", hrcm_vars[hrcm_bio]))
    # multiply by 10 to remove decimals and match global datasets
    tifX10<-og_tif*10
    # store bioclim variable name 
    new_name<-file_path_sans_ext(hrcm_vars[hrcm_bio])
    # save new raster as .tif file
    writeRaster(tifX10, filename = paste0(hrcm_path, new_name, ".tif"), 
                format = "GTiff", overwrite = TRUE)
  }
}

# run function and repeat for each spatial resolution of bioclims
edit_units_tif(c125)
edit_units_tif(c500)
edit_units_tif(f500)

### END HRCM UNITS ADJUSTMENTS ###