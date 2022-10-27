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
# remove Niihau
hi_map<-hi_map[which(hi_map$Island != "NI"),]
plot(hi_map)

# models
m_type<-c("global_notHI", "local_HI", "nested_HI")# c("global", "local", "nested")
m_type=paste0(m_type, "_models")

g_mod<-m_type[1]
l_mod<-m_type[2]
n_mod<-m_type[3]

# save outputs path
outDir<-paste0(wdDir, "combined_results/expert_analysis/")
dir.create(outDir, showWarnings = F, recursive = T)
dir.create(paste0(outDir, "suitability_maps/"), showWarnings = F)
dir.create(paste0(outDir, "binary_maps/"), showWarnings = F)

outDir2<-paste0(wdDir, "combined_results/combined_maps/")
dir.create(outDir2, showWarnings = F, recursive = T)
dir.create(paste0(outDir2, "suitability_maps/"), showWarnings = F)
dir.create(paste0(outDir2, "binary_maps/"), showWarnings = F)

av_sp<-c("Clid", "Falc", "Hedy", "Lant", "Leuc", "Meli", "More", 
         "Pani", "Pass", "Penc", "Pens", "Psid", "Seta", "Schi",  
         "Spha", "Mico", "Ulex")
# sp_nm %in% all_sp_nm
### SUMMARY PLOTS ###
par(mfrow = c(1, 1))

### MAPS ###
# # set color ramp (suitability and clipped)
# bin_col<-colorRampPalette(c("white", "forestgreen"))
# clip_col<-colorRampPalette(c("white", "dodgerblue3", "yellow", "red"))

# use vidris package for better color scale (colorblindness and gray scale)
# inferno() or magma()

# loop through and create species ranges comparisons per model 
eval_stat = eval_stats[1]
s=1
for (eval_stat in eval_stats){
  cat("doing ", eval_stat, "\n")
  for(s in 1:length(all_sp_nm)){  # set s = 1 for debugging
    # species name
    sp_nm=all_sp_nm[s]
    print(sp_nm)
    ### SUITABILITY PLOTS ###
    # clipped suitability to BIN range
    
    # vorsino raster
    v_clip<-raster(paste0(wdDir, "data/vorsino_files/critical_hab_analysis/baseline_habitat_analysis/Baseline_",
                          av_sp[s], "_Threshold_Scale_EM.tif"))
    # # set NA to 0 for unsuitable area
    # v_clip[is.na(v_clip)]<-0
    
    # global
    g_clip<-raster(paste0(wdDir, g_mod, "/output_rasters/", sp_nm, "_clipped_suitability_baseline_", eval_stat, "_wmean.tif"))
    g_clip[g_clip == 0]<-NA
    # local
    l_clip<-raster(paste0(wdDir, l_mod, "/output_rasters/", sp_nm, "_clipped_suitability_baseline_", eval_stat, "_wmean.tif"))
    l_clip[l_clip == 0]<-NA
    # nested
    n_clip<-raster(paste0(wdDir, n_mod, "/output_rasters/", sp_nm, "_clipped_suitability_baseline_", eval_stat, "_wmean.tif"))
    n_clip[n_clip == 0]<-NA
    
    # find max value
    zlim_clip<-max(c(summary(g_clip)[5], summary(l_clip)[5], summary(n_clip)[5]))
    #zlimMin_clip<-min(c(summary(g_clip)[1], summary(l_clip)[1], summary(n_clip)[1]))
    zlimMin_clip<-0 #min(c(summary(g_clip)[1], summary(l_clip)[1], summary(n_clip)[1]))
    
    # create image file
    png(paste0(outDir, "suitability_maps/", eval_stat, "_", sp_nm, "_baseline_clipped_suitability.png"), 
        width = 1500, height = 1200, units = "px", pointsize = 12, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # suitability model results
    plot(crop(g_clip, hi_map), main = "Model A", box = F, axes = F, legend = F,  
         #zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
         zlim = c(zlimMin_clip, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    plot(crop(l_clip, hi_map), main = "Model B", box = F, axes = F, legend = F,  
         zlim = c(zlimMin_clip, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    plot(crop(n_clip, hi_map), main = "Model C", box = F, axes = F, legend = F,  
         zlim = c(zlimMin_clip, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c("white", rev(inferno(4))), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # vorsino model
    plot(crop(v_clip, hi_map)/1000, main = "Model D", box = F, axes = F, legend = F,
         zlim = c(0, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    # add species name as main title 
    title(sub("_", " ", sp_nm), outer = T, cex.main = 2) #, line = 2)
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
    
    ######################
    #only global regional and nested
    # create image file
    png(paste0(outDir2, "suitability_maps/", eval_stat, "_", sp_nm, "_baseline_clipped_suitability.png"), 
        width = 1500, height = 600, units = "px", pointsize = 26, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(1, 3), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # suitability model results
    plot(crop(g_clip, hi_map), main = "Global", box = F, axes = F, legend = F,  
         #zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
         zlim = c(zlimMin_clip, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    plot(crop(l_clip, hi_map), main = "Regional", box = F, axes = F, legend = F,  
         zlim = c(zlimMin_clip, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    plot(crop(n_clip, hi_map), main = "Nested", box = F, axes = F, legend = F,  
         zlim = c(zlimMin_clip, zlim_clip), col = rev(inferno(30))); plot(hi_map, add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c("white", rev(inferno(4))), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
    # add species name as main title 
    # title(sub("_", " ", sp_nm), outer = T, cex.main = 2) #, line = 2)
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
    
    ###################
    #calculate some range metrics
    
    #calculate area
    pixel_area_sqkm=summary(area(l_bin))[3]
    g_bin_sqkm_area=freq(g_bin)[1,2]*pixel_area_sqkm
    l_bin_sqkm_area=freq(l_bin)[1,2]*pixel_area_sqkm
    n_bin_sqkm_area=freq(n_bin)[1,2]*pixel_area_sqkm
    area_DF_row=data.frame(eval_stat, sp_nm, g_bin_sqkm_area, l_bin_sqkm_area, n_bin_sqkm_area)
    
    #calculate niche overlap
    #clipped suitability
    library(ENMeval)
    g_clip0=g_clip
    g_clip0[is.na(g_clip0)]=0
    
    l_clip0=l_clip
    l_clip0[is.na(l_clip0)]=0
    
    n_clip0=n_clip
    n_clip0[is.na(n_clip0)]=0
    
    GL_clipSuit_D=calc.niche.overlap(stack(g_clip0, l_clip0), "D")[2,1]
    GL_clipSuit_I=calc.niche.overlap(stack(g_clip0, l_clip0), "I")[2,1]
    LN_clipSuit_D=calc.niche.overlap(stack(l_clip0, n_clip0), "D")[2,1]
    LN_clipSuit_I=calc.niche.overlap(stack(l_clip0, n_clip0), "I")[2,1]
    GN_clipSuit_D=calc.niche.overlap(stack(g_clip0, n_clip0), "D")[2,1]
    GN_clipSuit_I=calc.niche.overlap(stack(g_clip0, n_clip0), "I")[2,1]
    
    #binary
    g_bin0=g_bin
    g_bin0[is.na(g_bin0)]=0
    
    l_bin0=l_bin
    l_bin0[is.na(l_bin0)]=0
    
    n_bin0=n_bin
    n_bin0[is.na(n_bin0)]=0
    
    GL_bin_D=calc.niche.overlap(stack(g_bin0, l_bin0), "D")[2,1]
    GL_bin_I=calc.niche.overlap(stack(g_bin0, l_bin0), "I")[2,1]
    LN_bin_D=calc.niche.overlap(stack(l_bin0, n_bin0), "D")[2,1]
    LN_bin_I=calc.niche.overlap(stack(l_bin0, n_bin0), "I")[2,1]
    GN_bin_D=calc.niche.overlap(stack(g_bin0, n_bin0), "D")[2,1]
    GN_bin_I=calc.niche.overlap(stack(g_bin0, n_bin0), "I")[2,1]
    
    overlap_DF_row=data.frame(eval_stat, sp_nm, 
                              GL_clipSuit_D, LN_clipSuit_D, GN_clipSuit_D, 
                              GL_clipSuit_I, LN_clipSuit_I, GN_clipSuit_I,
                              GL_bin_D, LN_bin_D, GN_bin_D, 
                              GL_bin_I, LN_bin_I, GN_bin_I)
    
    if (s==1 & eval_stat==eval_stats[1]){
      area_DF=area_DF_row
      overlap_DF=overlap_DF_row
    }else{
      area_DF=rbind(area_DF, area_DF_row)
      overlap_DF=rbind(overlap_DF, overlap_DF_row)
    }
    
    # create image file
    png(paste0(outDir, "binary_maps/", eval_stat, "_", sp_nm, "_baseline_occurrence.png"), 
        width = 1500, height = 1200, units = "px", pointsize = 12, bg = "white", res = 82)
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
    title(sub("_", " ", sp_nm), outer = T, cex.main = 2) #, line = 2)
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
    
    # create image file
    png(paste0(outDir2, "binary_maps/", eval_stat, "_", sp_nm, "_baseline_occurrence.png"), 
        width = 1500, height = 600, units = "px", pointsize = 26, bg = "white", res = 82)
    # format plot area
    par(mfrow = c(1, 3), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
    # suitability model results
    plot(crop(g_bin, hi_map), main = "Global", box = F, 
         # axes = F, legend = F, col = bin_col(2)); plot(hi_map, add = T)
         axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
    plot(crop(l_bin, hi_map), main = "Regional", box = F, 
         axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
    plot(crop(n_bin, hi_map), main = "Nested", box = F, 
         axes = F, legend = F, col = rev(inferno(4))); plot(hi_map, add = T)
    # LEGEND OPTION 1: bottom left
    legend("bottomleft", bty = "n", fill = c("white", rev(inferno(4))[3]), border = "black", cex = 1.25,
           legend = c("Unsuitable", "Suitable"))
    # add species name as main title 
    title(sub("_", " ", sp_nm), outer = T, cex.main = 2) #, line = 2)
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
}
write.csv(area_DF, paste0(outDir2, "projection_area_DF.csv"), row.names = F)
write.csv(overlap_DF, paste0(outDir2, "niche_overlap_metric_DF.csv"), row.names = F)

