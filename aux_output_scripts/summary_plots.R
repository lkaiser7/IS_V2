# Phase 1 SDMs for Invasive Species
# additional figures for model runs
# global, regional, & weighted outputs

### SET UP ####

# working directory
wdDir<-rootDir
setwd(wdDir)

# load packages
library(rgdal)
library(raster)
graphics.off()
# load map data
hi_map<-readOGR(paste0(mapDir, "Main_Hawaiian_Islands_simple3.shp"))
# hi_map<-readOGR(paste0(wdDir, "data/map_data"), "state_coast") # SAME
# remove Niihau
hi_map<-hi_map[which(hi_map$Island != "NI"),]
plot(hi_map)

# models
m_type<-c("global_notHI", "regional_HI", "nested_HI")

g_mod<-paste0(m_type[1], "_models/")
l_mod<-paste0(m_type[2], "_models/")
n_mod<-paste0(m_type[3], "_models/")

# save outputs path
outDir<-paste0(wdDir, "combined_results/summary_results/")
dir.create(outDir, showWarnings = F, recursive = T)
dir.create(paste0(outDir, "suitability_maps/"), showWarnings = F)
dir.create(paste0(outDir, "clipped_suit_maps/"), showWarnings = F)
### SUMMARY PLOTS ###
# compare model outputs 
par(mfrow = c(1, 1))

### MAPS ###

# set color ramp (suitability and clipped)
bin_col<-colorRampPalette(c("white", "forestgreen"))
suit_col<-colorRampPalette(c("dodgerblue3", "yellow", "red"))
clip_col<-colorRampPalette(c("white", "dodgerblue3", "yellow", "red"))

# loop through and create species ranges comparisons per model 
eval_stat = eval_stats[1]
s=1
for (eval_stat in eval_stats){
  cat("doing ", eval_stat, "\n")
  dir.create(paste0(outDir, "bin_", eval_stat, "_maps/"), showWarnings = F)
  for(s in 1:length(all_sp_nm)){  # set s = 1 for debugging
    # species name
    current_sp_nm=replace_spp_names(all_sp_nm[s])
    print(current_sp_nm)
    
    # species regional points
    sp_pts<-read.csv(paste0(wdDir, "data/merged_data/hi_data/", all_sp_nm[s], ".csv"))
    
    ### BIN ROC PLOTS ###
    # compare ranges of models 
    
    # global
    g_bin<-raster(paste0(wdDir, g_mod, "output_rasters/", all_sp_nm[s], "_BIN_baseline_", eval_stat, "_wmean.tif"))
    # regional
    l_bin<-raster(paste0(wdDir, l_mod, "output_rasters/", all_sp_nm[s], "_BIN_baseline_", eval_stat, "_wmean.tif"))
    # nested
    n_bin<-raster(paste0(wdDir, n_mod, "output_rasters/", all_sp_nm[s], "_BIN_baseline_", eval_stat, "_wmean.tif"))
    
    # create image file
    png(paste0(outDir, "bin_", eval_stat, "_maps/", all_sp_nm[s], "_BIN_baseline_", eval_stat, "_wmean.png"), 
        width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 4, 0))
    # binary model results
    plot(g_bin, main = "Global Model", box = F, axes = F, legend = F); plot(hi_map, add = T)
    plot(l_bin, main = "Regional Model", box = F, axes = F, legend = F); plot(hi_map, add = T)
    plot(n_bin, main = "Nested Model", box = F, axes = F, legend = F); plot(hi_map, add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", legend = c("Unsuitable", "Suitable"),
           bty = "n", fill = bin_col(2), border = "black", cex = 1.5)
    # species point data
    plot(hi_map, main = "Regional Occurrences")
    points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
           pch = 21, col = "black", bg = "red") #sp_pts$dataColor)
    # legend("bottomleft", legend = c("Global Data", "National Data", "Regional Data"),
    #        pch = 19, col = c("blue", "yellow3", "green4"), bty = "n", cex = 1.5)
    # legend("bottomleft", legend = c("Global Data", "National Data", "Regional Data"),
    #        pch = 1, col = "black", bty = "n", cex = 1.5)
    # add species name as main title 
    title(current_sp_nm, outer = T, cex.main = 2) #, line = 2)
    dev.off()
    
    ### SUITABILITY PLOTS ###
    # compare suitability of species ranges
    
    # global
    g_suit<-raster(paste0(wdDir, g_mod, "output_rasters/", all_sp_nm[s], "_suitability_baseline_", eval_stat, "_wmean.tif"))
    g_suit[g_suit == 0]<-NA
    # regional
    l_suit<-raster(paste0(wdDir, l_mod, "output_rasters/", all_sp_nm[s], "_suitability_baseline_", eval_stat, "_wmean.tif"))
    l_suit[l_suit == 0]<-NA
    # nested
    n_suit<-raster(paste0(wdDir, n_mod, "output_rasters/", all_sp_nm[s], "_suitability_baseline_", eval_stat, "_wmean.tif"))
    n_suit[n_suit == 0]<-NA
    
    # find max value
    zlim_suit<-max(c(summary(g_suit)[5], summary(l_suit)[5], summary(n_suit)[5]))
    
    # create image file
    png(paste0(outDir, "suitability_maps/", all_sp_nm[s], "_suitability_baseline_", eval_stat, "_wmean.png"), 
        width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # suitability model results
    plot(g_suit, main ="Global Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, add = T)
    plot(l_suit, main = "Regional Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, add = T)
    plot(n_suit, main = "Nested Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # species point data
    plot(hi_map, main = "Regional Occurrences")
    points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
           pch = 21, col = "black", bg = "red")
    # add species name as main title 
    title(current_sp_nm, outer = T, cex.main = 2) #, line = 2)
    dev.off()
    
    ### CLIPPED PLOTS ###
    # clipped suitabilty to BIN range
    
    # global
    # g_clip<-g_suit*g_bin
    g_clip<-raster(paste0(wdDir, g_mod, "output_rasters/", all_sp_nm[s], "_clipped_suitability_baseline_", eval_stat, "_wmean.tif"))
    # g_0<-g_clip; g_0[g_0 != 0]<-NA
    # regional
    # l_clip<-l_suit*l_bin
    l_clip<-raster(paste0(wdDir, l_mod, "output_rasters/", all_sp_nm[s], "_clipped_suitability_baseline_", eval_stat, "_wmean.tif"))
    # nested
    # n_clip<-n_suit*n_bin
    n_clip<-raster(paste0(wdDir, n_mod, "output_rasters/", all_sp_nm[s], "_clipped_suitability_baseline_", eval_stat, "_wmean.tif"))
    
    # find max value
    zlim_clip<-max(c(summary(g_clip)[5], summary(l_clip)[5], summary(n_clip)[5]))
    
    # create image file
    png(paste0(outDir, "clipped_suit_maps/", all_sp_nm[s], "_clipped_suitability_baseline_", eval_stat, "_wmean.png"), 
        width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # suitability model results
    plot(g_clip, main = "Global Model", box = F, axes = F, legend = F,  
         zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
    # plot(g_0, add = T, col = "white", box = F, axes = F, legend = F); plot(hi_map, add = T)
    plot(l_clip, main = "Regional Model", box = F, axes = F, legend = F,  
         zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
    plot(n_clip, main = "Nested Model", box = F, axes = F, legend = F,  
         zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # species point data
    plot(hi_map, main = "Regional Occurrences")
    points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
           pch = 21, col = "black", bg = "red")
    # add species name as main title 
    title(current_sp_nm, outer = T, cex.main = 2) #, line = 2)
    dev.off()
    
    #merged plots
    ### MERGED BIN AND SUITABILITY PLOTS ###
    # compare suitability of species ranges
    
    # global
    # g_suit
    # l_suit
    # n_suit
    
    g_bin2=g_bin
    l_bin2=l_bin
    n_bin2=n_bin
    g_bin2[g_bin2==0]=NA
    l_bin2[l_bin2==0]=NA
    n_bin2[n_bin2==0]=NA
    g_bin3=aggregate(g_bin2, 4, modal)
    l_bin3=aggregate(l_bin2, 4, modal)
    n_bin3=aggregate(n_bin2, 4, modal)
    g_bin_pol=rasterToPolygons(g_bin3, dissolve=T) 
    l_bin_pol=rasterToPolygons(l_bin3, dissolve=T)
    n_bin_pol=rasterToPolygons(n_bin3, dissolve=T)
    library(rgeos)
    #sum(sapply(g_bin_pol@polygons, function(y) nrow(y@Polygons[[2]]@coords))) #total vertices
    g_bin_pol_simple=gSimplify(g_bin_pol, 0.0005, topologyPreserve=TRUE)
    l_bin_pol_simple=gSimplify(l_bin_pol, 0.0005, topologyPreserve=TRUE)
    n_bin_pol_simple=gSimplify(n_bin_pol, 0.0005, topologyPreserve=TRUE)
    #sum(sapply(g_bin_pol_simple@polygons, function(y) nrow(y@Polygons[[2]]@coords))) #total vertices
    plot(g_bin_pol)
    plot(g_bin_pol_simple)
    plot(l_bin_pol_simple)
    plot(n_bin_pol_simple)
    
    ##########################
    #save for data release
    DRdir=paste0(outDir, "combo_bin_and_suitability_maps_3_panel/", "data_release_files/")
    dir.create(DRdir, showWarnings = F, recursive = T) #all_sp_nm[s]
    writeRaster(x = g_suit, filename = paste0(DRdir, current_sp_nm, " global_suitability.tif"), overwrite=T)#, gdal=c("COMPRESS=LZW"))
    writeRaster(g_bin2, filename = paste0(DRdir, current_sp_nm, " global_binary_range.tif"), overwrite=T)#, gdal=c("compress=lzw"))
    terra::writeVector(terra::vect(g_bin_pol_simple), filename = paste0(DRdir, current_sp_nm, " global_binary_range.gpkg"), overwrite=T)#, overwrite=T)

    writeRaster(l_suit, filename = paste0(DRdir, current_sp_nm, " regional_suitability.tif"), overwrite=T)#, gdal=c("compress=lzw"))
    writeRaster(l_bin2, filename = paste0(DRdir, current_sp_nm, " regional_binary_range.tif"), overwrite=T)#, gdal=c("compress=lzw"))
    terra::writeVector(terra::vect(l_bin_pol_simple), filename = paste0(DRdir, current_sp_nm, " regional_binary_range.gpkg"), overwrite=T)
    
    writeRaster(n_suit, filename = paste0(DRdir, current_sp_nm, " nested_suitability.tif"), overwrite=T)#, gdal=c("compress=lzw"))
    writeRaster(n_bin2, filename = paste0(DRdir, current_sp_nm, " nested_binary_range.tif"), overwrite=T)#, gdal=c("compress=lzw"))
    terra::writeVector(terra::vect(n_bin_pol_simple), filename = paste0(DRdir, current_sp_nm, " nested_binary_range.gpkg"), overwrite=T)
    
    # find max value
    zlim_suit<-max(c(summary(g_suit)[5], summary(l_suit)[5], summary(n_suit)[5]))
    
    dir.create(paste0(outDir, "combo_bin_and_suitability_maps/"), showWarnings=F)
    # create image file
    png(paste0(outDir, "combo_bin_and_suitability_maps/", all_sp_nm[s], "_combo_baseline_", eval_stat, "_wmean.png"), 
        width = 870*2, height = 607*2, units = "px", pointsize = 24, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    library(scales)
    # suitability model results
    plot(g_suit, main ="Global Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(g_bin_pol_simple, col = alpha("black", 0.25), add = T)
    plot(l_suit, main = "Regional Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(l_bin_pol_simple, col = alpha("black", 0.25), add = T)
    plot(n_suit, main = "Nested Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(n_bin_pol_simple, col = alpha("black", 0.25), add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # species point data
    plot(hi_map, main = "Regional Occurrences")
    points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
           pch = 21, col = "black", bg = "red")
    # add species name as main title 
    title(current_sp_nm, outer = T, cex.main = 2) #, line = 2)
    dev.off()
    
    ########################################
    #modified 3 panel combo plot
    dir.create(paste0(outDir, "combo_bin_and_suitability_maps_3_panel_mod/"), showWarnings=F)
    # create image file
    png(paste0(outDir, "combo_bin_and_suitability_maps_3_panel_mod/", all_sp_nm[s], "_combo_baseline_", eval_stat, "_wmean.png"), 
        width = 870*1, height = 607*3, units = "px", pointsize = 24, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(3, 1), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    library(scales)
    # suitability model results
    plot(g_suit, main ="Global Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(g_bin_pol_simple, col = alpha("black", 0.25), add = T)
    plot(l_suit, main = "Regional Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(l_bin_pol_simple, col = alpha("black", 0.25), add = T)
    plot(n_suit, main = "Nested Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(n_bin_pol_simple, col = alpha("black", 0.25), add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # # species point data
    # plot(hi_map, main = "Regional Occurrences")
    # points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
    #        pch = 21, col = "black", bg = "red")
    # add species name as main title 
    title(current_sp_nm, outer = T, cex.main = 2) #, line = 2)
    dev.off()
    
    ######################
    #3 panel combo plot
    dir.create(paste0(outDir, "combo_bin_and_suitability_maps_3_panel/"), showWarnings=F)
    # create image file
    png(paste0(outDir, "combo_bin_and_suitability_maps_3_panel/", all_sp_nm[s], "_combo_baseline_", eval_stat, "_wmean.png"), 
        width = 870*2, height = 500, units = "px", pointsize = 24, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(1, 3), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # suitability model results
    plot(g_suit, main ="Global Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(g_bin_pol_simple, col = alpha("black", 0.25), add = T)
    plot(l_suit, main = "Regional Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(l_bin_pol_simple, col = alpha("black", 0.25), add = T)
    plot(n_suit, main = "Nested Model", box = F, axes = F, legend = F, 
         zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, lwd=2, border="grey", add = T); plot(n_bin_pol_simple, col = alpha("black", 0.25), add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # # species point data
    # plot(hi_map, main = "Regional Occurrences")
    # points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
    #        pch = 21, col = "black", bg = "red")
    # add species name as main title 
    title(current_sp_nm, outer = T, cex.main = 2) #, line = 2)
    dev.off()
    
  }
}

# reset plot
graphics.off()
#dev.new()

