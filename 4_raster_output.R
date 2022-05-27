### BIOMOD2: creates rasters ###
### see source script inputs ###
### step 4:

##################
##### SET UP #####
##################
setwd(rootDir)
# load necessary packages 
library(biomod2)
library(stringr)
library(colorRamps)
library(raster)
######################
##### FUNCTIONS ######
######################

# function to save display and save raster image
save_raster_fx = function(raster_img, out_nm) {
  # save name of jpeg file to be created
  jpeg_name = paste0(out_nm, ".jpg")
  # create blank jpeg file 
  jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
       width = 10, height = 8, quality = 90, bg = "white")
  # plot image
  plot(raster_img)
  # save image file
  dev.off() 
  
  # save name of tiff file to be created
  tiff_name = paste0(out_nm, ".tif")
  # save name of raster file to be created
  out_raster_name = paste0(out_nm, ".tif")
  # create blank tiff file 
  tiff(tiff_name, res = 300, units = "in", pointsize = 12, 
       width = 10, height = 8, compression = "lzw")
  # plot image
  plot(raster_img)
  # save image file
  dev.off()  
  
  # save raster data file
  writeRaster(raster_img, out_raster_name, format="GTiff", overwrite = TRUE)        
}


# function
Process_raster_data_BadtoGood = function(raster_var, out_nm, min_lim = NULL, 
                                         max_lim = NULL, mask_data = NULL){
  # save name of jpeg file
  jpeg_name = paste0(out_nm, ".jpg")
  # save name of raster file to be created
  out_raster_name = paste0(out_nm, ".tif")
  # create blank jpeg file
  jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
       width = 14, height = 14, quality = 90, bg = "white")
  
  # add minimum if default is set to NULL
  if (is.null(min_lim)){
    # find minimum from raster
    min_lim = minValue(raster_var)
  }
  # add maximum if default is set to NULL
  if (is.null(max_lim)){
    # find maximum from raster
    max_lim = maxValue(raster_var)
  }
  if (maxValue(raster_var)>max_lim) max_lim=maxValue(raster_var)
  
  # set color ramp palette
  col5<-colorRampPalette(c('red', 'gray96', 'darkgreen')) 
  # plot raster
  plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
       axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
       legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
       axis.args = list(at = seq(min_lim, max_lim, (max_lim-min_lim)/10),
                      labels = seq(min_lim, max_lim, (max_lim-min_lim)/10)))
  
  # check if mask data is NULL 
  if (!is.null(mask_data)){
    # add mask to plot
    plot(mask_data, add = TRUE)    
  }
  # save jpeg image file
  dev.off()  
  # save raster data file
  writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite = TRUE)
}

# function
Process_raster_data_NeutraltoGood = function(raster_var, out_nm, min_lim = NULL, 
                                             max_lim = NULL, mask_data = NULL){
  # save name of jpeg file
  jpeg_name = paste0(out_nm, ".jpg")
  # save name of raster file to be created
  out_raster_name = paste0(out_nm, ".tif")
  # create blank jpeg file
  jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
       width = 14, height = 14, quality = 90, bg = "white")
  
  # add minimum if default is set to NULL
  if (is.null(min_lim)){
    # find minimum from raster
    min_lim = minValue(raster_var)
  }
  # add maximum if default is set to NULL
  if (is.null(max_lim)){
    # find maximum from raster
    max_lim = maxValue(raster_var)
  }
  if (maxValue(raster_var)>max_lim) max_lim=maxValue(raster_var)
  
  # set color ramp palette
  col5<-colorRampPalette(c('gray96', 'darkgreen')) 
  # plot raster
  plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
       axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
       legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
       axis.args = list(at = seq(min_lim, max_lim, (max_lim-min_lim)/10),
                      labels = seq(min_lim, max_lim, (max_lim-min_lim)/10)))
  
  # check if mask data is NULL 
  if (!is.null(mask_data)){
    # add mask to plot
    plot(mask_data, add = TRUE)    
  }
  # save jpeg image file
  dev.off()  
  # save raster data file
  writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite=TRUE)
}

# function
Process_raster_data_NeutraltoBad = function(raster_var, out_nm, min_lim = NULL, 
                                            max_lim = NULL, mask_data = NULL){
  # save name of jpeg file
  jpeg_name = paste0(out_nm, ".jpg")
  # save name of raster file to be created
  out_raster_name = paste0(out_nm, ".tif")
  # create blank tiff file
  # create blank jpeg file
  jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
       width = 14, height = 14, quality = 90, bg = "white")
  
  # add minimum if default is set to NULL
  if (is.null(min_lim)){
    # find minimum from raster
    min_lim = minValue(raster_var)
  }
  # add maximum if default is set to NULL
  if (is.null(max_lim)){
    # find maximum from raster
    max_lim = maxValue(raster_var)
  }

  # set color ramp palette
  col5<-colorRampPalette(c('gray96', 'red'))
  # plot raster
  plot(raster_var, col = col5(n = 99), breaks = seq(min_lim, max_lim, length.out = 100), 
       axes = FALSE, box = FALSE, legend = TRUE, legend.width = 1, legend.shrink = 0.75,
       legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
       axis.args = list(at = seq(min_lim, max_lim, (max_lim-min_lim)/10),
                      labels = seq(min_lim, max_lim, (max_lim-min_lim)/10)))
  
  # check if mask data is NULL 
  if (!is.null(mask_data)){
    # add mask to plot
    plot(mask_data, add = TRUE)    
  }
  # save jpeg image file
  dev.off()  
  # save raster data file
  writeRaster(raster_var, out_raster_name, format = "GTiff", overwrite = TRUE)
}

# set mask layer to shapefile of the main Hawaiian islands
mask_layer = shapefile(paste0(mapDir, "Main_Hawaiian_Islands_simple3.shp"))

#########################
##### RUN FUNCTIONS #####
#########################

# create directory for rasters
dir.create(paste0(project_path, 'output_rasters/'), showWarnings = FALSE)
# create directory for main raster result outputs
dir.create(paste0(project_path, 'output_rasters/main/'), showWarnings = FALSE)

# set species name to first species
sp_nm = all_sp_nm[1]

# select first evaluation statistic
eval_stat = spp_ensemble_eval_stats[3]
# loop through all evaluation statistics
for (eval_stat in spp_ensemble_eval_stats){
  # loop through all species 
  for (sp_nm in all_sp_nm){
    # convert species name to character
    sp_nm = as.character(sp_nm)  
    # store species name as character
    sp_nm0 = sp_nm
    # replace species naming convention of "_" with "."
    sp_nm = str_replace_all(sp_nm, "_", ".")
    
    # print sign posting of ongoing modeling per species 
    cat('\n',sp_nm,'raster output creation')
    
    # store file output name for species response zones per evaluation and ensemble type
    out_nm = paste0(project_run, '/output_rasters/main/', sp_nm0, "_response_zones_", 
                    eval_stat, "_", spp_ensemble_type, "_", comp_projects[2])
    # store file output name as tiff
    out_raster_name00 = paste0(out_nm,".tif")
    
    # if species file already exists, skip to end, otherwise begin raster creation
    if (file.exists(out_raster_name00) == FALSE | overwrite == 1){
      
      #######################
      ##### BINARY MAPS #####
      #######################
      
      # store raster names
      raster_names = c("EM_suitability1", "EM_suitability2")
      # store bin names
      raster_names_bin = c("EM_BIN1", "EM_BIN2")
      
      # reset i counter to 1
      i = 1
      # loop through 2
      for (i in 1:projections_to_run){
        # store individual project name
        proj_nm = comp_projects[i]
        # store individual raster name
        raster_name = raster_names[i]
        # store individual bin name
        raster_name_bin = raster_names_bin[i]
        
        # store first file name for raster .grd file
        file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                            proj_nm, "_", sp_nm, "_ensemble.tif")
        # store temporary raster stack 
        temp_raster = stack(file_name1)
        #names(temp_raster)
        # find integer of band of total consensus weighted mean ROC from ensemble models 
        band_n = which(names(temp_raster) == paste0(sp_nm, '_EM', spp_ensemble_type,
                                                    'By', eval_stat, 
                                                    '_mergedAlgo_mergedRun_mergedData'))        
        # assign stored raster name to selected band with weighted mean ROC values
        assign(raster_name, raster(temp_raster, layer = band_n)/1000)

        # store first bin name for raster .grd file
        file_name1_bin = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                                proj_nm, "_", sp_nm, "_ensemble_", eval_stat, "bin.tif")
        # store temporary raster stack 
        temp_raster_bin = stack(file_name1_bin)  
        
        # find integer of desired band
        band_n = which(names(temp_raster) == paste0(sp_nm,'_EM',spp_ensemble_type,
                                                    'By',eval_stat, 
                                                    '_mergedAlgo_mergedRun_mergedData')) 
        # assign stored raster bin name to selected band 
        assign(raster_name_bin, raster(temp_raster_bin, layer = band_n))
        
        # run if true in source script
        if (plot_spp_ensemble_CV){
          # find integer of desired band
          band_n = which(names(temp_raster) == paste0(sp_nm,'_EMcvBy',eval_stat,
                                                      '_mergedAlgo_mergedRun_mergedData'))
          # assign stored raster name to selected band
          assign(paste0(raster_name, "_CV"), raster(temp_raster, layer = band_n)/1000)
          
          # create output suitability rasters for CV 
          out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_", "suitability_CV_", 
                          proj_nm, "_", eval_stat, "_", spp_ensemble_type)
          # run function from above
          Process_raster_data_NeutraltoBad(get(paste0(raster_name, "_CV")),
                                           out_nm, mask_data = mask_layer)          
        }
        
        # create output suitability rasters for each evaluation statistic
        out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_", "suitability_", 
                        proj_nm, "_", eval_stat, "_", spp_ensemble_type)
        # run function from above
        Process_raster_data_NeutraltoGood(get(raster_name), mask_data = mask_layer,
                                          out_nm, min_lim = 0, max_lim = 1)
        
        # create output bin rasters for each evaluation statistic   
        out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_", "BIN_", 
                        proj_nm, "_", eval_stat, "_", spp_ensemble_type)
        # run function from above to save raster image
        save_raster_fx(get(raster_name_bin), out_nm)
      }
      
      # print sign posting of loading of rasters complete
      cat('\n done with loading baseline and future rasters for', sp_nm)
      
      # masked species ensemble map suitability 1
      masked_suitability1 = EM_BIN1*EM_suitability1
      # create clipped suitability map for each evaluation statistic
      out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_", "clipped_suitability_", 
                      comp_projects[1], "_", eval_stat, "_", spp_ensemble_type)
      # run function from above to save raster image
      # save_raster_fx(masked_suitability1, out_nm)
      # run function from above
      Process_raster_data_NeutraltoGood(masked_suitability1, out_nm,min_lim = 0, 
                                        max_lim = 1, mask_data = mask_layer)
      
      # skip following if only creating baseline rasters 
      if (projections_to_run == 2){
        # masked species ensemble map suitability 2
        masked_suitability2 = EM_BIN2*EM_suitability2
        # create clipped suitability map for each evaluation statistic
        out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_", "clipped_suitability_", 
                        comp_projects[2], "_", eval_stat, "_", spp_ensemble_type)
        # run function from above to save raster image
        # save_raster_fx(masked_suitability2, out_nm)
        # run function from above
        Process_raster_data_NeutraltoGood(masked_suitability2, out_nm, min_lim = 0, 
                                          max_lim = 1, mask_data = mask_layer)
        
        ### FIX ###
        # change in masked species ensemble suitability maps - DIFFERENT RESOLUTIONS???
        suitability_change = EM_suitability2-EM_suitability1
        # create change in suitability map
        out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_", "suitability_change_", 
                        eval_stat, "_", spp_ensemble_type)
        # run function from above to save raster image
        # save_raster_fx(suitability_change, out_nm)
        # run function from above
        Process_raster_data_BadtoGood(suitability_change, out_nm, min_lim = -1, 
                                      max_lim = 1, mask_data = mask_layer)
        
        # create temporary vector for ensemble bins
        jnk = EM_BIN2*10
        # calculate sum of ensemble bins
        BIN_dif = EM_BIN1+jnk
        # create vector 
        m = c(9.9, 10.1, 3, 10.9, 11.1, 2)
        # create matrix
        rclmat = matrix(m,  ncol = 3,  byrow = TRUE)
        # reclassify response zones based on bins
        resp_zone = reclassify(BIN_dif,  rclmat)
        
        # create color palette
        mypalette_numbers = c(0, 1, 2, 3)
        # select colors for palettes
        mypalette = c("Grey", "Red", "Green", "Yellow")
        # create vector of response zone names
        resp_zone_names0 = c("Lost", "Overlap", "Gained")
        
        # run if true in source script
        if (masked_spp_ensemble_map){
          # set current baseline mask
          current_mask = EM_suitability1 > minValue(EM_suitability1) 
          # store future analog climates raster
          analog_cc_loc = paste0(sp_nm0, "_analog_climates2100.tif")
          # load analog climates raster
          analog_cc = raster(analog_cc_loc)
          
          # combine into one mask
          all_mask = analog_cc*2 + current_mask*4 
          #1 cur, 2 ang, 4 hab, 3 cur/ang, 6 ang/hab, 7 cur/ang/hab
          # calculate cumulative mask
          multi_mask = current_mask*analog_cc
          # create mask of response zones
          masked_resp_zone = resp_zone*multi_mask
          
          # print sign posting of created masks for species
          cat('\n created mask for', sp_nm)
          
          # save name of jpeg file for mask 
          jpeg_name = paste0(project_run, '/output_rasters/main/', sp_nm0, "_mask.jpg")
          # create blank jpeg file
          jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
               width = 10, height = 8, quality = 90, bg = "white")
          # plot mask
          plot(all_mask)
          # save tiff file
          dev.off()
          
          # masked bin comparison rasters
          out_nm = paste0(project_run, '/output_rasters/main/', sp_nm0, "_response_zones_masked_",
                          eval_stat, "_", spp_ensemble_type, "_", comp_projects[2])
          # save name of jpeg file
          jpeg_name = paste0(out_nm, ".jpg")
          # save name of raster file to be created
          out_raster_name = paste0(out_nm, ".tif")
          
          # create temporary object of response zones from mask
          jnk = unique(masked_resp_zone)
          # create palette of colors
          graph_palette = mypalette_numbers
          # select minimum range for present zone
          zones_present = jnk[jnk > 0]
          # select ,aximum range for present zone
          zones_present = zones_present[zones_present <= 3]
          # assign palette colors to zones
          resp_zone_colors = mypalette[zones_present + 1]
          # assign names to zones
          resp_zone_names = resp_zone_names0[zones_present]
          # select unique zone numbers for palette from temporary object
          mypalette_numbers_selected = mypalette[jnk + 1]
          
          # create blank jpeg file
          jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
               width = 10, height = 8, quality = 90, bg = "white")
          # plot masked response zone raster
          plot(masked_resp_zone, col = mypalette_numbers_selected, legend = FALSE)
          # add legend
          legend("bottomleft", legend = resp_zone_names, col = resp_zone_colors, pch = 16)
          # save raster image file
          dev.off()  
          # save raster data file
          writeRaster(masked_resp_zone, out_raster_name, format = "GTiff", overwrite = TRUE)
          
          # combine with cumulative mask
          future_bin_with_mask = multi_mask*EM_BIN2
          # add suitability mask
          future_suitability_with_mask = multi_mask*EM_suitability2
          
          # create output suitability rasters for each image
          out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_suitability_future_masked_", 
                          eval_stat, "_", spp_ensemble_type, "_", comp_projects[2])
          # save name of jpeg file
          jpeg_name = paste0(out_nm, ".jpg")
          # save name of raster file to be created
          out_raster_name = paste0(out_nm, ".tif")
          # create blank jpeg file
          jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
               width = 10, height = 8, quality = 90, bg = "white")
          # plot future suitability mask
          plot(future_suitability_with_mask)
          # save raster image file
          dev.off()
          # save raster data file        
          writeRaster(future_suitability_with_mask, out_raster_name, 
                      format = "GTiff", overwrite = TRUE)
          
          # masked binary maps per evaluation statistic
          out_nm = paste0(project_run, '/output_rasters/', sp_nm0, "_BIN_future_masked_", 
                          eval_stat, "_", spp_ensemble_type, "_", comp_projects[2])
          # save name of jpeg file
          jpeg_name = paste0(out_nm, ".jpg")
          # save name of raster file to be created
          out_raster_name = paste0(out_nm, ".tif")
          # create blank jpeg file
          jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
               width = 10, height = 8, quality = 90, bg = "white")
          # plot binary map with mask
          plot(future_bin_with_mask)
          # save raster image file
          dev.off()
          # save raster data file
          writeRaster(future_bin_with_mask, out_raster_name,  
                      format = "GTiff", overwrite = TRUE)
        }
        
        
        # response zones comarison rasters per evaluation statistic
        out_nm = paste0(project_run, '/output_rasters/main/', sp_nm0, "_response_zones_",
                        eval_stat, "_", spp_ensemble_type, "_", comp_projects[2])
        # save name of jpeg file
        jpeg_name = paste0(out_nm, ".jpg")
        # save name of raster file to be created
        out_raster_name00 = paste0(out_nm, ".tif")
        
        # create temporary object of response zones
        jnk = unique(resp_zone)
        # create palette of colors
        graph_palette = mypalette_numbers
        # select minimum range for present zone
        zones_present = jnk[jnk > 0]
        # select ,aximum range for present zone
        zones_present = zones_present[zones_present <= 3]
        # assign palette colors to zones
        resp_zone_colors = mypalette[zones_present + 1]
        # assign names to zones
        resp_zone_names = resp_zone_names0[zones_present]
        # select unique zone numbers for palette from temporary object
        mypalette_numbers_selected = mypalette[jnk + 1]    
        
        # create blank jpeg file
        jpeg(jpeg_name, res = 300, units = "in", pointsize = 12, 
             width = 10, height = 8, quality = 90, bg = "white")
        # plot response zones raster
        plot(resp_zone,  main = sp_nm0, col = mypalette_numbers_selected, legend = FALSE)
        # add legend
        legend("bottomleft", legend = resp_zone_names, col = resp_zone_colors, pch = 16)
        # save raster image file
        dev.off()
        # save raster data file 
        writeRaster(resp_zone, out_raster_name00, format = "GTiff", overwrite = TRUE)
      }
    
    # otherwise if raster outputs are already done
    }else{
      # then print sign posting of completed per species
      cat('\n', sp_nm, 'already calculated')
    }
  }
}

##############################
##### END RASTER OUTPUTS #####
##############################