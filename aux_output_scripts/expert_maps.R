# Phase 1 SDMs for Invasive Species
# additional figure for expert analysis
# global, local, nested & original maps

### SET UP ####

# working directory
wdDir<-rootDir
setwd(wdDir)

# load packages
library(rgdal)
library(raster)
library(viridis)

# load map data
hi_map<-readOGR(paste0(mapDir, "Main_Hawaiian_Islands_simple3.shp"))
# hi_map<-readOGR(paste0(wdDir, "data/map_data"), "state_coast") # SAME
# remove Niihau
hi_map<-hi_map[which(hi_map$Island != "NI"),]
plot(hi_map)

# models
m_type<-c("global_notHI", "local_HI", "nested_HI")# c("global", "local", "nested")
# # set model type
# mod<-m_type[3]
# select main or LRG species
# sdm species
sp_nm<-c("Clidemia_hirta", "Falcataria_moluccana", "Hedychium_gardnerianum", 
         "Lantana_camara", "Leucaena_leucocephala", "Melinis_minutiflora",
         "Morella_faya", "Panicum_maximum", "Passiflora_tarminiana", 
         "Pennisetum_clandestinum", "Pennisetum_setaceum", "Psidium_cattleianum", 
         "Schinus_terebinthifolius", "Setaria_palmifolia", 
         "Cyathea_cooperi", "Miconia_calvescens", "Ulex_europaeus")
av_sp<-c("Clid", "Falc", "Hedy", "Lant", "Leuc", "Meli", "More", 
         "Pani", "Pass", "Penc", "Pens", "Psid", "Schi", "Seta", 
         "Spha", "Mico", "Ulex")

# Cyathea.cooperi run for both in global and local models but use LRG outputs

g_mod<-m_type[1]
l_mod<-m_type[2]
n_mod<-m_type[3]


# save outputs path
outDir<-paste0(wdDir, "expert_analysis/")

### SUMMARY PLOTS ###
# compare model outputs 
par(mfrow = c(1, 1))

### MAPS ###

# # set color ramp (suitability and clipped)
# bin_col<-colorRampPalette(c("white", "forestgreen"))
# clip_col<-colorRampPalette(c("white", "dodgerblue3", "yellow", "red"))

# use vidris package for better color scale (colorblindness and gray scale)
# inferno() or magma()

# loop through and create species ranges comparisons per model 
for(s in 1:length(sp_nm)){  # set s = 1 for debugging
  # species name
  print(sp_nm[s])
  print(av_sp[s])
  
  ### SUITABILITY PLOTS ###
  # clipped suitability to BIN range
  
  # vorsino raster
  v_clip<-raster(paste0(wdDir, "data/vorsino_files/Crit_hab_analysis/Baseline_",
                        av_sp[s], "_Threshold_Scale_EM.tif"))
  # # set NA to 0 for unsuitable area
  # v_clip[is.na(v_clip)]<-0
  
  # global
  g_clip<-raster(paste0(wdDir, g_mod, "output_rasters/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.tif"))
  g_clip[g_clip == 0]<-NA
  # local
  l_clip<-raster(paste0(wdDir, l_mod, "output_rasters/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.tif"))
  l_clip[l_clip == 0]<-NA
  # nested
  n_clip<-raster(paste0(wdDir, n_mod, "output_rasters/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.tif"))
  n_clip[n_clip == 0]<-NA
  
  # find max value
  zlim_clip<-max(c(summary(g_clip)[5], summary(l_clip)[5], summary(n_clip)[5]))
  
  # create image file
  png(paste0(outDir, "suitability_maps/", sp_nm[s], "_baseline_clipped_suitability.png"), 
      width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
  # format plot area
  par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
  # suitability model results
  plot(crop(g_clip, hi_map), main = "Model A", box = F, axes = F, legend = F,  
       #zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
       zlim = c(0, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
  plot(crop(l_clip, hi_map), main = "Model B", box = F, axes = F, legend = F,  
       zlim = c(0, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
  plot(crop(n_clip, hi_map), main = "Model C", box = F, axes = F, legend = F,  
       zlim = c(0, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
  # LEGEND OPTION 1: bottom left
  legend("bottomleft", bty = "n", fill = c("white", rev(inferno(4))), border = "black", cex = 1.25,
         legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
  # vorsino model
  plot(crop(v_clip, hi_map)/1000, main = "Model D", box = F, axes = F, legend = F,
       zlim = c(0, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
  # add species name as main title 
  title(sub("_", " ", sp_nm[s]), outer = T, cex.main = 2) #, line = 2)
  # # LEGEND OPTION 2: under title
  # # overlay plot area to add legend 
  # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 1.5, 0), mar = c(0, 0, 0, 0), new = TRUE)
  # plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  # # add legend
  # # strwidth(c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
  # legend("top", c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"),
  #        text.width = c(0.25, 0.3, 0.35, 0.45, 0.25),
  #        xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", fill = clip_col(5), border = "black", cex = 1.25)
  # save plot
  dev.off()
  
  ### BINARY PLOTS ###
  # range of occurrences
  
  # vorsino raster
  v_bin<-v_clip
  v_bin[v_bin > 0]<-1
  
  # global
  g_bin<-g_clip
  g_bin[g_bin > 0]<-1
  # local
  l_bin<-l_clip
  l_bin[l_bin > 0]<-1
  # nested
  n_bin<-n_clip
  n_bin[n_bin > 0]<-1
  
  # create image file
  png(paste0(outDir, "binary_maps/", sp_nm[s], "_baseline_occurrance.png"), 
      width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
  # format plot area
  par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
  # suitability model results
  plot(crop(g_bin, hi_map), main = "Model A", box = F, 
       # axes = F, legend = F, col = bin_col(2)); plot(hi_map, add = T)
       axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
  plot(crop(l_bin, hi_map), main = "Model B", box = F, 
       axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
  plot(crop(n_bin, hi_map), main = "Model C", box = F, 
       axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
  # LEGEND OPTION 1: bottom left
  legend("bottomleft", bty = "n", fill = c("white", rev(inferno(4))[3]), border = "black", cex = 1.25,
         legend = c("Unsuitable", "Suitable"))
  # vorsino model
  plot(crop(v_bin, hi_map), main = "Model D", box = F, 
       axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
  # add species name as main title 
  title(sub("_", " ", sp_nm[s]), outer = T, cex.main = 2) #, line = 2)
  # # LEGEND OPTION 2: under title
  # # overlay plot area to add legend 
  # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 1.5, 0), mar = c(0, 0, 0, 0), new = TRUE)
  # plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  # # add legend
  # # strwidth(c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
  # legend("top", c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"),
  #        text.width = c(0.25, 0.3, 0.35, 0.45, 0.25),
  #        xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", fill = clip_col(5), border = "black", cex = 1.25)
  # save plot
  dev.off()
  
}
