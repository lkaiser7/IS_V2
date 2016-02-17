### BIOMOD2: mapped ensemble ###
### see source script inputs ###
### step 7:

##################
##### SET UP #####
##################

# load necessary packages
library(biomod2)
library(stringr)

# set working directory to project path
setwd(project_path)

# create directory for tables in project folder 
dir.create("tables/", showWarnings = FALSE)

# set climate data directory to baseline bioclims
clim_data_dir = biobaseRun 

# store vector of variable to include in map
vars = c("Eval stat", "Species", "baseline_area", "future_area", "% change", 
         "lost", "kept", "gained", 
         "baseline_suitability", "future_suitability", "% suitability change")

####################################################
##### FUNCTIONS TO MAP SPECIES ENSEMBLE RANGES #####
####################################################

# load bioclim raster to build mask
island_mask = raster(paste0(clim_data_dir, "bio1.tif"))
# create binary mask from bioclim raster
island_mask[!is.na(island_mask)] = 1

# function to return raster with larger spatial extent that has max and min coordinates
extend_raster = function(raster_lyr){
  # set any missing or NA values to 0
  raster_lyr[raster_lyr == NA] = 0
  # extend raster layer by entire spatial extent
  raster_lyr = extend(raster_lyr, island_mask, value = 0)
  # combine extended extent with binary raster mask
  expanded_raster_lyr = raster_lyr*island_mask
  # return the expanded and combined raster layer
  return(expanded_raster_lyr)  
}

# FUNCTION
Process_raster_data_BadtoGood = 
  function(raster_var, out_nm, min_lim = NULL, max_lim = NULL, mask_data = NULL){
    
    # store file output name for jpeg image
    jpeg_name = paste0(out_nm, ".jpg")
    # store file output name for tiff file
    out_raster_name = paste0(out_nm, ".tif")
    # create blank jpeg file
    jpeg(jpeg_name, width = 14, height = 14, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    
    # set minimum limit to raster
    if (is.null(min_lim)){
      min_lim = minValue(raster_var)
    }
    # set maximum limit to raster
    if (is.null(max_lim)){
      max_lim=maxValue(raster_var)
    }
  
    # load colorRamps package from library
    library(colorRamps)
    # set color ramp to desired palette
    col5<-colorRampPalette(c('red', 'gray96', 'darkgreen'))  
  
    # plot raster
    plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
         axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
         legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
         axis.args = list(at = seq(min_lim, max_lim, (max_lim-min_lim)/10),
                          labels = seq(min_lim, max_lim, (max_lim-min_lim)/10)))
  
    # add mask to map
    if (!is.null(mask_data)){
      plot(mask_data, add = TRUE)  
    }
    
    # save jpeg image file
    dev.off()
    # save raster as tiff file
    writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite = TRUE)
}

# FUNCTION
Process_raster_data_NeutraltoGood = 
  function(raster_var, out_nm, min_lim = NULL, max_lim = NULL, mask_data = NULL){
    
    jpeg_name = paste0(out_nm, ".jpg")
    out_raster_name = paste0(out_nm, ".tif")
    jpeg(jpeg_name, width = 14, height = 14, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    
    if (is.null(min_lim)){
      min_lim = minValue(raster_var)
    }
    if (is.null(max_lim)){
      max_lim = maxValue(raster_var)
    }
    
    library(colorRamps)
    col5<-colorRampPalette(c('gray96', 'darkgreen'))
    
    plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
         axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
         legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
         axis.args = list(at = seq(min_lim,max_lim, (max_lim-min_lim)/10),
                        labels = seq(min_lim,max_lim, (max_lim-min_lim)/10)))
    
    if (!is.null(mask_data)){
      plot(mask_data, add = TRUE)    
    }
    
    dev.off()  
    writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite = TRUE)
}

# FUNCTION
Process_raster_data_NeutraltoGood_W_overlay = 
  function(raster_var, out_nm, min_lim = NULL, max_lim = NULL, 
           mask_data = NULL, overlay_data = NULL){
    
    jpeg_name=paste(out_nm, ".jpg", sep = "")
    out_raster_name=paste(out_nm, ".tif", sep = "")
    jpeg(jpeg_name, width = 14, height = 14, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
  
    if (is.null(min_lim)){
      min_lim = minValue(raster_var)
    }
    if (is.null(max_lim)){
      max_lim = maxValue(raster_var)
    }
  
    library(colorRamps)
    col5<-colorRampPalette(c('gray96', 'darkgreen'))
 
    plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
         axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
         legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
         axis.args = list(at = seq(min_lim, max_lim, (max_lim-min_lim)/10),
                        labels = seq(min_lim, max_lim, (max_lim-min_lim)/10)))
  
    if (!is.null(mask_data)){
      plot(mask_data,add = TRUE)    
    }
    if (!is.null(overlay_data)){
      plot(overlay_data, add = T, border = "red", lwd = 3)    
    }
   
    dev.off()  
    writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite = TRUE)
}

# FUNCTION
Process_raster_data_NeutraltoBad = function(raster_var, out_nm, min_lim = NULL, 
                                            max_lim = NULL, mask_data = NULL){
  
  jpeg_name = paste0(out_nm, ".jpg")
  out_raster_name = paste0(out_nm, ".tif")
  jpeg(jpeg_name, width = 14, height = 14, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  
  if (is.null(min_lim)){
    min_lim = minValue(raster_var)
  }
  if (is.null(max_lim)){
    max_lim = maxValue(raster_var)
  }
  
  library(colorRamps)
  col5<-colorRampPalette(c('gray96', 'red'))
  
  plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
       axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
       legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
       axis.args = list(at = seq(min_lim, max_lim, (max_lim-min_lim)/10),
                      labels = seq(min_lim, max_lim, (max_lim-min_lim)/10)))
  
  if (!is.null(mask_data)){
    plot(mask_data, add = TRUE)    
  }
  
  dev.off()  
  writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite = TRUE)
}

######################################################
##### RUN FUNCTIONS TO MAP SPECIES RANGE OUTPUTS #####
######################################################

# set species name to first species
sp_nm = all_sp_nm[1]
# set ensemble evaluation statistic to first stat
eval_stat = spp_ensemble_eval_stats[1]

# definte output raster folder
out_ras<-paste0(project_path, "output_rasters/")

# load mask layer shapefile
mask_layer = shapefile(paste0(mapDir, "Main_Hawaiian_Islands_simple3.shp"))

# loop through all ensemble evaluation statistics 
for (eval_stat in spp_ensemble_eval_stats){
  # current suitability, future suitability
  # masked current suitability, masked future suitability
  # current bin, future bin
  # suitability_delta (change)
  # lost_range, kept_range, gained_range
  
  # store raster mask as temporary file
  jnk = island_mask
  # set mask values to 0
  jnk[jnk == 1] = 0
  
  # store desired variable outcomes as blank masks
  spp_em_masked_current_suitability = jnk
  spp_em_masked_future_suitability = jnk
  spp_em_masked_current_suitability_CV = jnk
  spp_em_masked_future_suitability_CV = jnk
  spp_em_current_suitability = jnk
  spp_em_current_suitability_bin = jnk
  spp_em_future_suitability = jnk
  spp_em_future_suitability_bin = jnk
  spp_em_current_bin = jnk
  spp_em_future_bin = jnk
  spp_em_suitability_delta = jnk
  spp_em_lost_range = jnk
  spp_em_kept_range = jnk
  spp_em_gained_range = jnk

  # loop through all species
  for (sp_nm in all_sp_nm){
    # convert species name to character
    sp_nm = as.character(sp_nm)  
    # store species name for file names
    sp_nm0 = sp_nm
    # replace species naming convention of "_" with "." 
    sp_nm = str_replace_all(sp_nm,"_", ".")
    
    # print sign posting of ongoing mapping per species
    cat('\n',sp_nm0,'ensemble mapping...')
    
    # load raster masks for baseline, future, and suitability 
    response_raster = raster(paste0(out_ras, "main/", sp_nm0, "_response_zones_", eval_stat,
                                    "_", spp_ensemble_type, "_", comp_projects[2], ".tif"))  
    baseline_masked_suitability = raster(paste0(out_ras, sp_nm0, "_", "clipped_suitability_",
                                                "baseline", "_", eval_stat, "_", spp_ensemble_type, ".tif"))
    future_masked_suitability = raster(paste0(out_ras, sp_nm0, "_", "clipped_suitability_",
                                              comp_projects[2], "_", eval_stat, "_", spp_ensemble_type, ".tif"))
    # load ensemble suitability rasters
    baseline_suitability = raster(paste0(out_ras, sp_nm0, "_", "suitability_", "baseline", "_", 
                                         eval_stat, "_", spp_ensemble_type, ".tif"))
    baseline_suitability_bin = baseline_suitability > 0
    future_suitability = raster(paste0(out_ras, sp_nm0, "_", "suitability_", comp_projects[2],
                                       "_", eval_stat, "_", spp_ensemble_type, ".tif"))
    future_suitability_bin = future_suitability > 0
    # load cv suitability rasters
    baseline_suitability_CV = raster(paste0(out_ras, sp_nm0, "_", "suitability_CV_", "baseline",
                                            "_", eval_stat, "_", spp_ensemble_type, ".tif"))
    future_suitability_CV = raster(paste0(out_ras, sp_nm0, "_", "suitability_CV_", 
                                          comp_projects[2], "_", eval_stat, "_", spp_ensemble_type, ".tif"))
    # load bin rasters and suitability change
    current_bin = raster(paste0(out_ras, sp_nm0, "_", "BIN_", "baseline", "_", eval_stat, "_",
                                spp_ensemble_type, ".tif"))
    future_bin = raster(paste0(out_ras, sp_nm0, "_", "BIN_", comp_projects[2], "_", eval_stat,
                               "_", spp_ensemble_type, ".tif"))
    suitability_delta = raster(paste0(out_ras, sp_nm0, "_", "suitability_change_", eval_stat,
                                      "_", spp_ensemble_type, ".tif"))
    # set range values for lost, kept, and gained ranges
    lost_range = response_raster == 1
    kept_range = response_raster == 2
    gained_range = response_raster == 3
    
    # align and extend raster masks to account for all water
    baseline_masked_suitability = extend_raster(baseline_masked_suitability)
    future_masked_suitability = extend_raster(future_masked_suitability)  
    baseline_suitability = extend_raster(baseline_suitability)
    baseline_suitability_bin = extend_raster(baseline_suitability_bin)
    future_suitability = extend_raster(future_suitability)
    future_suitability_bin = extend_raster(future_suitability_bin)
    baseline_suitability_CV = extend_raster(baseline_suitability_CV)
    future_suitability_CV = extend_raster(future_suitability_CV)
    current_bin = extend_raster(current_bin)
    future_bin = extend_raster(future_bin)
    suitability_delta = extend_raster(suitability_delta)
    lost_range = extend_raster(lost_range)
    kept_range = extend_raster(kept_range)
    gained_range = extend_raster(gained_range)
    
    # raster calculations
    spp_em_masked_current_suitability = spp_em_masked_current_suitability + baseline_masked_suitability
    spp_em_masked_future_suitability = spp_em_masked_future_suitability + future_masked_suitability
    spp_em_current_suitability = spp_em_current_suitability + baseline_suitability
    spp_em_current_suitability_bin = spp_em_current_suitability_bin + baseline_suitability_bin
    spp_em_future_suitability = spp_em_future_suitability + future_suitability
    spp_em_future_suitability_bin = spp_em_future_suitability_bin + future_suitability_bin
    spp_em_masked_current_suitability_CV = spp_em_masked_current_suitability_CV + baseline_suitability_CV
    spp_em_masked_future_suitability_CV = spp_em_masked_future_suitability_CV + future_suitability_CV
    spp_em_current_bin = spp_em_current_bin + current_bin
    spp_em_future_bin = spp_em_future_bin + future_bin
    spp_em_suitability_delta = spp_em_suitability_delta + suitability_delta
    spp_em_lost_range = spp_em_lost_range + lost_range
    spp_em_kept_range = spp_em_kept_range + kept_range
    spp_em_gained_range = spp_em_gained_range + gained_range
  }
  
  # create folder to save ensemble maps
  dir.create('output_rasters/spp_ensembles/', showWarnings = FALSE)
  # run functions with loaded rasters and save outputs
  Process_raster_data_NeutraltoGood(spp_em_current_suitability/spp_em_current_suitability_bin, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_current_suitability_', eval_stat), 
                                    max_lim = 1, min_lim = 0, mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_future_suitability/spp_em_future_suitability_bin, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_', comp_projects[2], '_suitability_', eval_stat), 
                                    max_lim = 1, min_lim = 0, mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_masked_current_suitability/spp_em_current_bin, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_masked_current_suitability_', eval_stat), 
                                    max_lim = 1, min_lim = 0, mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_masked_future_suitability/spp_em_future_bin, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_masked_', comp_projects[2], '_suitability_', eval_stat),
                                    max_lim = 1, min_lim = 0, mask_data = mask_layer)
  Process_raster_data_NeutraltoBad(spp_em_masked_current_suitability_CV, 
                                   paste0(out_ras, 'spp_ensembles/spp_em_current_suitability_CV_', eval_stat), 
                                   mask_data = mask_layer)
  Process_raster_data_NeutraltoBad(spp_em_masked_future_suitability_CV, 
                                   paste0(out_ras, 'spp_ensembles/spp_em_', comp_projects[2], '_suitability_CV_', eval_stat), 
                                   mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_current_bin, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_current_bin_', eval_stat), 
                                    mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_future_bin, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_', comp_projects[2], '_bin_', eval_stat), 
                                    mask_data = mask_layer)
  
  # calculate average species ensemble suitability
  avg_spp_em_suitability_delta = (spp_em_future_suitability/spp_em_future_suitability_bin) - (spp_em_current_suitability/spp_em_current_suitability_bin)
  
  # set max value
  max_val = max(abs(c(cellStats(avg_spp_em_suitability_delta, max), cellStats(avg_spp_em_suitability_delta, min))))
  # run functions with loaded rasters
  Process_raster_data_BadtoGood(avg_spp_em_suitability_delta, 
                                paste0(out_ras, 'spp_ensembles/spp_em_suitability_avg_delta_', comp_projects[2],'_', eval_stat), 
                                max_lim = max_val, min_lim = -max_val, mask_data = mask_layer)
  
  # reset max value
  max_val=max(abs(c(cellStats(spp_em_suitability_delta, max), cellStats(spp_em_suitability_delta, min))))
  # run functions with loaded rasters
  Process_raster_data_BadtoGood(spp_em_suitability_delta, 
                                paste0(out_ras, 'spp_ensembles/spp_em_suitability_delta_', comp_projects[2], '_', eval_stat), 
                                max_lim = max_val, min_lim = -max_val, mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_lost_range, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_lost_range_', comp_projects[2], '_', eval_stat), 
                                    mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_kept_range, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_kept_range_', comp_projects[2], '_', eval_stat), 
                                    mask_data = mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_gained_range, 
                                    paste0(out_ras, 'spp_ensembles/spp_em_gained_range_', comp_projects[2], '_', eval_stat), 
                                    mask_data = mask_layer)
  
  # shapefile figures with results for publication
  prot_areas = shapefile(paste0("Y:/PICCC_analysis/FB_analysis/habitat_analysis/",
                                "protected_areas_20100331_simpleWGS1984.shp"))
  # run function with loaded rasters with overlay of protected areas
  Process_raster_data_NeutraltoGood_W_overlay(spp_em_kept_range, 
                                              paste0(out_ras, 'spp_ensembles/spp_em_kept_ranges_w_prot_areas_', comp_projects[2], '_', eval_stat), 
                                              mask_data = mask_layer, overlay_data = prot_areas)
}

################################
##### END ENSEMBLE MAPPING #####
################################