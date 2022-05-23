# Phase 1 SDMs for Invasive Species
# additional figures for model runs
# global, local, & weighted outputs

### SET UP ####

# working directory
wdDir<-rootDir
setwd(wdDir)

# load packages
library(rgdal)
library(raster)

# load map data
hi_map<-readOGR(paste0(mapDir, "Main_Hawaiian_Islands_simple3.shp"))
# hi_map<-readOGR(paste0(wdDir, "data/map_data"), "state_coast") # SAME
# remove Niihau
hi_map<-hi_map[which(hi_map$Island != "NI"),]
plot(hi_map)

# models
m_type<-c("global_notHI", "local_HI", "nested_HI")
# # set model type
# mod<-m_type[3]
# select main or LRG species
LRG = F

# sdm species
sp_sm<-c("Clidemia_hirta", "Falcataria_moluccana", "Hedychium_gardnerianum", 
         "Lantana_camara", "Leucaena_leucocephala", "Melinis_minutiflora",
         "Morella_faya", "Panicum_maximum", "Passiflora_tarminiana", 
         "Pennisetum_clandestinum", "Pennisetum_setaceum", "Psidium_cattleianum", 
         "Schinus_terebinthifolius", "Setaria_palmifolia")

# LRG species
sp_lrg<-c("Cyathea_cooperi", "Miconia_calvescens", "Ulex_europaeus")
sp_lrg<-c("Cyathea_cooperi", "Ulex_europaeus")
# Cyathea.cooperi run for both in global and local models but use LRG outputs

# # select model inputs
# if(mod == "global" & LRG){
#   # model path
#   m_path<-paste0(mod, "_notHI_models_LRG/")
#   # species names
#   sp_nm<-sp_lrg
# }else if(mod == "global" & !LRG){
#   # model path
#   m_path<-paste0(mod, "_notHI_models/")
#   # species names
#   sp_nm<-sp_sm
# }else if(mod != "global" & LRG){
#   # model path
#   m_path<-paste0(mod, "_HI_models_LRG/")
#   # species names
#   sp_nm<-sp_lrg
# }else{
#   # model path
#   m_path<-paste0(mod, "_HI_models/")
#   # species names
#   sp_nm<-sp_sm
# }
# 
# # check path and species
# m_path; length(sp_nm)
# 
# # output folder paths
# rastDir<-paste0(wdDir, m_path, "output_rasters/")
# plotDir<-paste0(wdDir, m_path, "outputs/")

# if(LRG){
#   # model paths
#   g_mod<-paste0(m_type[1], "_notHI_models_LRG/")
#   l_mod<-paste0(m_type[2], "_HI_models_LRG/")
#   n_mod<-paste0(m_type[3], "_HI_models_LRG/")
#   # species names
#   sp_nm<-sp_lrg
#   
# }else{
  # model paths
  g_mod<-paste0(m_type[1], "/")
  l_mod<-paste0(m_type[2], "/")
  n_mod<-paste0(m_type[3], "/")
  # species names
  sp_nm<-sp_sm
#   
# }

# save outputs path
outDir<-paste0(wdDir, "summary_results/")

### SUMMARY PLOTS ###
# compare model outputs 
par(mfrow = c(1, 1))

### MAPS ###

# set color ramp (suitability and clipped)
bin_col<-colorRampPalette(c("white", "forestgreen"))
suit_col<-colorRampPalette(c("dodgerblue3", "yellow", "red"))
clip_col<-colorRampPalette(c("white", "dodgerblue3", "yellow", "red"))

# loop through and create species ranges comparisons per model 
for(s in 1:length(sp_nm)){  # set s = 1 for debugging
  # species name
  print(sp_nm[s])
  
  # species local points
  sp_pts<-read.csv(paste0(wdDir, "data/hi_data/", sp_nm[s], ".csv"))
  # # add color column
  # colIndex<-c("G", "B", "L")
  # colValues<-c("blue", "yellow3", "green4")
  # sp_pts$dataColor<-colValues[match(sp_pts$dataset, colIndex)]
  # # head(sp_pts)
  
  ### BIN ROC PLOTS ###
  # compare ranges of models 
  
  # global
  g_bin<-raster(paste0(wdDir, g_mod, "output_rasters/", sp_nm[s], "_BIN_baseline_ROC_wmean.tif"))
  # local
  l_bin<-raster(paste0(wdDir, l_mod, "output_rasters/", sp_nm[s], "_BIN_baseline_ROC_wmean.tif"))
  # nested
  n_bin<-raster(paste0(wdDir, n_mod, "output_rasters/", sp_nm[s], "_BIN_baseline_ROC_wmean.tif"))
  
  # create image file
  png(paste0(outDir, "bin_ROC_maps/", sp_nm[s], "_BIN_baseline_ROC_wmean.png"), 
      width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
  # format plot area
  par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
  # par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 4, 0))
  # binary model results
  plot(g_bin, main = "Global Model", box = F, axes = F, legend = F); plot(hi_map, add = T)
  plot(l_bin, main = "Local Model", box = F, axes = F, legend = F); plot(hi_map, add = T)
  plot(n_bin, main = "Nested Model", box = F, axes = F, legend = F); plot(hi_map, add = T)
  # LEGEND OPTION 1: bottom left
  legend("bottomleft", legend = c("Unsuitable", "Suitable"),
         bty = "n", fill = bin_col(2), border = "black", cex = 1.5)
  # species point data
  plot(hi_map, main = "Local Occurrences")
  points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
         pch = 21, col = "black", bg = "red") #sp_pts$dataColor)
  # legend("bottomleft", legend = c("Global Data", "National Data", "Local Data"),
  #        pch = 19, col = c("blue", "yellow3", "green4"), bty = "n", cex = 1.5)
  # legend("bottomleft", legend = c("Global Data", "National Data", "Local Data"),
  #        pch = 1, col = "black", bty = "n", cex = 1.5)
  # add species name as main title 
  title(sub("_", " ", sp_nm[s]), outer = T, cex.main = 2) #, line = 2)
  # # LEGEND OPTION 2: under title
  # # overlay plot area to add legend 
  # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 1.5, 0), mar = c(0, 0, 0, 0), new = TRUE)
  # plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  # # add legend
  # legend("top", c("Unsuitable", "Suitable"), xpd = TRUE, horiz = TRUE, 
  #        inset = c(0, 0), bty = "n", fill = bin_col(2), border = "black", cex = 1.25)
  # save plot
  dev.off()
  
  ### SUITABILITY PLOTS ###
  # compare suitability of species ranges
  
  # global
  g_suit<-raster(paste0(wdDir, g_mod, "output_rasters/", sp_nm[s], "_suitability_baseline_ROC_wmean.tif"))
  g_suit[g_suit == 0]<-NA
  # local
  l_suit<-raster(paste0(wdDir, l_mod, "output_rasters/", sp_nm[s], "_suitability_baseline_ROC_wmean.tif"))
  l_suit[l_suit == 0]<-NA
  # nested
  n_suit<-raster(paste0(wdDir, n_mod, "output_rasters/", sp_nm[s], "_suitability_baseline_ROC_wmean.tif"))
  n_suit[n_suit == 0]<-NA
  
  # find max value
  zlim_suit<-max(c(summary(g_suit)[5], summary(l_suit)[5], summary(n_suit)[5]))
  
  # create image file
  png(paste0(outDir, "suitability_maps/", sp_nm[s], "_suitability_baseline_ROC_wmean.png"), 
      width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
  # format plot area
  par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
  # suitability model results
  plot(g_suit, main ="Global Model", box = F, axes = F, legend = F, 
       zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, add = T)
  plot(l_suit, main = "Local Model", box = F, axes = F, legend = F, 
       zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, add = T)
  plot(n_suit, main = "Nested Model", box = F, axes = F, legend = F, 
       zlim = c(0, zlim_suit), col = suit_col(30)); plot(hi_map, add = T)
  # LEGEND OPTION 1: bottom left
  legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
         legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
  # species point data
  plot(hi_map, main = "Local Occurrences")
  points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
         pch = 21, col = "black", bg = "red")
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
  
  ### CLIPPED PLOTS ###
  # clipped suitabilty to BIN range
  
  # global
  # g_clip<-g_suit*g_bin
  g_clip<-raster(paste0(wdDir, g_mod, "output_rasters/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.tif"))
  # g_0<-g_clip; g_0[g_0 != 0]<-NA
  # local
  # l_clip<-l_suit*l_bin
  l_clip<-raster(paste0(wdDir, l_mod, "output_rasters/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.tif"))
  # nested
  # n_clip<-n_suit*n_bin
  n_clip<-raster(paste0(wdDir, n_mod, "output_rasters/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.tif"))
  
  # find max value
  zlim_clip<-max(c(summary(g_clip)[5], summary(l_clip)[5], summary(n_clip)[5]))
  
  # create image file
  png(paste0(outDir, "clipped_suit_maps/", sp_nm[s], "_clipped_suitability_baseline_ROC_wmean.png"), 
      width = 870, height = 607, units = "px", pointsize = 12, bg = "white", res = 82)
  # format plot area
  par(mfrow = c(2, 2), mar = c(0.75, 0, 0.75, 0), oma = c(0, 0, 2, 0))
  # suitability model results
  plot(g_clip, main = "Global Model", box = F, axes = F, legend = F,  
       zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
  # plot(g_0, add = T, col = "white", box = F, axes = F, legend = F); plot(hi_map, add = T)
  plot(l_clip, main = "Local Model", box = F, axes = F, legend = F,  
       zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
  plot(n_clip, main = "Nested Model", box = F, axes = F, legend = F,  
       zlim = c(0, zlim_clip), col = clip_col(30)); plot(hi_map, add = T)
  # LEGEND OPTION 1: bottom left
  legend("bottomleft", bty = "n", fill = c(clip_col(5)), border = "black", cex = 1.25,
         legend = c("Unsuitable", "Low Suitability", "Moderate Suitability", "Suitable", "High Suitability"))
  # species point data
  plot(hi_map, main = "Local Occurrences")
  points(sp_pts$decimalLongitude, sp_pts$decimalLatitude, 
         pch = 21, col = "black", bg = "red")
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

# reset plot
dev.off()
dev.new()

### CURVES ###

# library (ggplot2)

# load response curve data
# global
g_roc<-read.csv(paste0(wdDir, g_mod, "outputs/all_eval_mat_ROC.csv"))
# local
l_roc<-read.csv(paste0(wdDir, l_mod, "outputs/all_eval_mat_ROC.csv"))
# nested
n_roc<-read.csv(paste0(wdDir, n_mod, "outputs/all_eval_mat_ROC.csv"))

head(n_roc)

# ### EVAL METRIC COMPARISON ###
# # function to calculate mean
# calculate_eval_mean_eval_metric = function(metric_df, model_eval_metric = "MAXENT.Phillips"){
#   metric_df_maxent_mean = metric_df[metric_df$rownames.Spp_eval. == model_eval_metric, ]
#   metric_df_maxent_mean = metric_df_maxent_mean[,-2]
#   jnk = apply(metric_df_maxent_mean[,-1], 1, mean, na.rm = T)
#   mean_eval_metric = data.frame(species = metric_df_maxent_mean[, 1], eval_metric = jnk)
#   return(mean_eval_metric)
# }
# 
# # MAXENT
# mean_global_maxent_eval = calculate_eval_mean_eval_metric(g_roc)
# mean_local_maxent_eval = calculate_eval_mean_eval_metric(l_roc)
# mean_nest_maxent_eval = calculate_eval_mean_eval_metric(n_roc)
# # GBM
# mean_global_GBM_eval = calculate_eval_mean_eval_metric(g_roc, model_eval_metric="GBM")
# mean_local_GBM_eval = calculate_eval_mean_eval_metric(l_roc, model_eval_metric="GBM")
# mean_nest_GBM_eval = calculate_eval_mean_eval_metric(n_roc, model_eval_metric="GBM")
# 
# # BY SPECIES
# mean_maxent_eval_df = merge(mean_global_maxent_eval, mean_local_maxent_eval, by = "species")
# mean_GBM_eval_df = merge(mean_global_GBM_eval, mean_local_GBM_eval, by = "species")
# # add nested data
# if(LRG){
#   mean_maxent_eval_df<-merge(mean_maxent_eval_df, mean_nest_maxent_eval)
#   mean_GBM_eval_df<-merge(mean_GBM_eval_df, mean_nest_GBM_eval)
# }else{
#   mean_maxent_eval_df<-merge(mean_maxent_eval_df[-nrow(mean_maxent_eval_df),], mean_nest_maxent_eval, by = "species")
#   mean_GBM_eval_df<-merge(mean_GBM_eval_df[-nrow(mean_GBM_eval_df),], mean_nest_GBM_eval, by = "species")
# }
# # rename function output table
# names(mean_maxent_eval_df) = c("Species", "Global", "Local", "Nested")
# names(mean_GBM_eval_df) = c("Species", "Global", "Local", "Nested")
# 
# # global v local
# plot(mean_maxent_eval_df$Global, mean_maxent_eval_df$Local)
# plot(mean_GBM_eval_df$Global, mean_GBM_eval_df$Local)
# # global v nested
# plot(mean_maxent_eval_df$Global, mean_maxent_eval_df$Nested)
# plot(mean_GBM_eval_df$Global, mean_GBM_eval_df$Nested)
# # local v nested
# plot(mean_maxent_eval_df$Local, mean_maxent_eval_df$Nested)
# plot(mean_GBM_eval_df$Local, mean_GBM_eval_df$Nested)

# loop through and compare model runs per species 
for(r in 1:length(sp_nm)){  # set r = 1 for debugging

  ### RESPONSE CURVE PLOTS ####
  # compare niches of species 
  par(mfrow = c(1, 1))

  ### MEAN - single values only ###
  # # global
  # g_roc<-read.csv(paste0(wdDir, g_mod, "outputs/all_eval_mean_mat_ROC.csv"))
  # # local
  # l_roc<-read.csv(paste0(wdDir, l_mod, "outputs/all_eval_mean_mat_ROC.csv"))
  # # nested
  # n_roc<-read.csv(paste0(wdDir, n_mod, "outputs/all_eval_mean_mat_ROC.csv"))
  
  # loop ticker
  print(sp_nm[r])
  
  # select species records
  g_sp<-cbind(model = "global", g_roc[which(g_roc[1] == sp_nm[r]),])
  l_sp<-cbind(model = "local", l_roc[which(l_roc[1] == sp_nm[r]),])
  n_sp<-cbind(model = "nested", n_roc[which(n_roc[1] == sp_nm[r]),])
  
  # combine records to loop
  all_sp<-rbind(g_sp, l_sp, n_sp)
  # all_sp
  # sp_means<-cbind(all_sp[,1:3], rowMeans(all_sp[,-1:-3]))
  
  # # all records
  # maxent g_sp[1, -1:-2]
  # gbm g_sp[2, -1:-2]
  
  # # line plots
  # plot(as.numeric(g_sp[1, -1:-2]), type = 'o')
  # lines(as.numeric(g_sp[2, -1:-2]), type = 'o', col = 'red')
  
  # # density plots
  # plot(density(as.numeric(g_sp[1, -1:-2])))
  # plot(density(as.numeric(g_sp[2, -1:-2])))
  
  # create image file
  png(paste0(outDir, "all_roc_curves/", sp_nm[r], "_all_eval_mat_ROC.png"),
      width = 680, height = 880, units = "px", pointsize = 18)
  # set plot 
  par(mfrow = c(3, 1))
  # loop through each model run
  for(m in 1:3){  # set m = 1 for debugging 
    # select model type
    m_name<-m_type[m]
    
    # select species by model type
    m_sp<-all_sp[which(all_sp$model == m_name),]
    
    # set up data variables 
    x_max<-as.numeric(m_sp[1, -1:-3])
    max_mean<-mean(x_max)
    x_gbm<-as.numeric(m_sp[2, -1:-3])
    gbm_mean<-mean(x_gbm)
    
    # # check number of bins
    # sqrt(length(x_max))
    # sqrt(length(x_gbm))
    
    h_max<-hist(x_max, breaks = 10, plot = F)
    xfit_max<-seq(min(x_max), max(x_max), length = 40)
    yfit_max<-dnorm(xfit_max, mean = mean(x_max), sd = sd(x_max))
    yfit_max<-yfit_max*diff(h_max$mids[1:2])*length(x_max)
    # lines(xfit_max, yfit_max, col = "blue", lwd = 2)
    
    h_gbm<-hist(x_gbm, breaks = 10, plot = F)  
    xfit_gbm<-seq(min(x_gbm), max(x_gbm), length = 40)
    yfit_gbm<-dnorm(xfit_gbm, mean = mean(x_gbm), sd = sd(x_gbm))
    yfit_gbm<-yfit_gbm*diff(h_gbm$mids[1:2])*length(x_gbm)
    # lines(xfit_gbm, yfit_gbm, col = "blue", lwd = 2)
    
    # set lims
    d_xlim<-c(min(xfit_max, xfit_gbm), max(xfit_max, xfit_gbm))
    d_ylim<-c(min(yfit_max, yfit_gbm), max(yfit_max, yfit_gbm))
    
    # normal curve plots
    plot(xfit_max, yfit_max, xlim = d_xlim, ylim = d_ylim,
         col = "blue", lwd = max_mean*4, type = 'l', #lwd = 1.5, type = 'o', 
         xlab = "Distribution", ylab = "Frequency",
         main = paste(sub("_", " ", sp_nm[r]), m_name, "ROC Model Evaluation Metrics"))
    lines(xfit_gbm, yfit_gbm, col = "red", lwd = gbm_mean*4, type = 'l') #lwd = 1.5, type = 'o')
    # legend("topleft", legend = c("MaxEnt", "GBM"), 
    #        col = c("blue", "red"), pch = 1, bty = "n")
    legend("bottomright", legend = c("MaxEnt", "GBM"), horiz = TRUE, 
           fill = c("blue", "red"), bty = "n", inset = c(0, 0.99), xpd = TRUE)
    legend("top", text.font = 4, bty = "n", inset = c(0, 1.12), xpd = T, 
           paste("MaxEnt Mean =", round(max_mean, 3), "GBM Mean =", round(gbm_mean, 3)))
  }
  # save plot
  dev.off()
  rm(m)

  # NOT ENOUGH DATA POINTS TO PLOT DISTRIBUTION CURVES
  # # full runs only
  # full_nm<-grep("Full", names(g_sp), value = T)
  # # maxent g_sp[1 , full_nm]
  # # gbm g_sp[2 , full_nm]
  # 
  # # # line plots
  # # plot(as.numeric(g_sp[1 , full_nm]), type = 'o', ylim = c(0.96, 0.98))
  # # lines(as.numeric(g_sp[2 , full_nm]), type = 'o', ylim = c(0.96, 0.98), col = 'red')
  # 
  # # # density plots
  # # plot(density(as.numeric(g_sp[1 , full_nm])), xlim = c(0.965, 0.985), ylim = c(0,700))
  # # lines(density(as.numeric(g_sp[2 , full_nm])), col = 'red')
  # 
  # # create image file
  # png(paste0(outDir, "full_roc_curves/", sp_nm[r], "_full_eval_mat_ROC.png"), 
  #     width = 680, height = 880, units = "px", pointsize = 18)
  # # set plot 
  # par(mfrow = c(3, 1))
  # # loop through each model run
  # for(m in 1:3){  # set m = 1 for debugging 
  #   # select model type
  #   m_name<-m_type[m]
  #   
  #   # select species by model type
  #   m_sp<-all_sp[which(all_sp$model == m_name),]
  #   
  #   # set up data variables - full runs only
  #   x_max<-as.numeric(m_sp[1, full_nm])
  #   x_gbm<-as.numeric(m_sp[2, full_nm])
  #   
  #   # # check number of bins
  #   # sqrt(length(x_max))
  #   # sqrt(length(x_gbm))
  #   
  #   h_max<-hist(x_max, breaks = 3, plot = F)
  #   xfit_max<-seq(min(x_max), max(x_max), length = 40)
  #   yfit_max<-dnorm(xfit_max, mean = mean(x_max), sd = sd(x_max))
  #   yfit_max<-yfit_max*diff(h_max$mids[1:2])*length(x_max)
  #   # lines(xfit_max, yfit_max, col = "blue", lwd = 2)
  #   
  #   h_gbm<-hist(x_gbm, breaks = 3, plot = F)  
  #   xfit_gbm<-seq(min(x_gbm), max(x_gbm), length = 40)
  #   yfit_gbm<-dnorm(xfit_gbm, mean = mean(x_gbm), sd = sd(x_gbm))
  #   yfit_gbm<-yfit_gbm*diff(h_gbm$mids[1:2])*length(x_gbm)
  #   # lines(xfit_gbm, yfit_gbm, col = "blue", lwd = 2)
  #   
  #   # set lims
  #   d_xlim<-c(min(xfit_max, xfit_gbm), max(xfit_max, xfit_gbm))
  #   d_ylim<-c(min(yfit_max, yfit_gbm), max(yfit_max, yfit_gbm))
  #   
  #   # normal curve plots
  #   plot(xfit_max, yfit_max, xlim = d_xlim, ylim = d_ylim,
  #        col = "blue", lwd = 1.5, type = 'o', 
  #        main = paste(sp_nm[r], m_name, "ROC Model Evaluation Metrics"))
  #   lines(xfit_gbm, yfit_gbm, col = "red", lwd = 1.5, type = 'o')
  #   legend("topleft", legend = c("MaxEnt", "GBM"), 
  #          col = c("blue", "red"), pch = 1, bty = "n")
  # }
  # # save plot
  # dev.off()
  
  ### VARIABLE IMPORTANCE CURVE PLOTS ####
  # compare bioclims per model for species 
  par(mfrow = c(1, 1))
  
  # global
  g_var<-read.csv(paste0(wdDir, g_mod, "outputs/all_VariImp.csv"))
  # local
  l_var<-read.csv(paste0(wdDir, l_mod, "outputs/all_VariImp.csv"))
  # nested
  n_var<-read.csv(paste0(wdDir, n_mod, "outputs/all_VariImp.csv"))
  
  ### MEAN - single values only ###
  # # global
  # g_var<-read.csv(paste0(wdDir, g_mod, "outputs/all_VariImp_mean.csv"))
  # # local
  # l_var<-read.csv(paste0(wdDir, l_mod, "outputs/all_VariImp_mean.csv"))
  # # nested
  # n_var<-read.csv(paste0(wdDir, n_mod, "outputs/all_VariImp_mean.csv"))
  
  # head(g_var)
  
  # select species
  g_sp<-g_var[which(g_var[1] == sp_nm[r]),]
  l_sp<-l_var[which(l_var[1] == sp_nm[r]),]
  n_sp<-n_var[which(n_var[1] == sp_nm[r]),]
  
  # combine records to loop
  all_sp<-list(g_sp, l_sp, n_sp)
  # all_sp
  
  # list variables 
  var_nm<-g_sp[, 2]
  
  # list column names by model type
  max_names<-grep("MAXENT", names(g_sp), value = T)
  gbm_names<-grep("GBM", names(g_sp), value = T)
  
  # max_vals<-g_sp[,max_names]
  # gbm_vals<-g_sp[,gbm_names]
  
  # # 1, 7, 12, 15
  # g_mean<-rowMeans(g_sp[,max_names])

  # # line plots
  # plot(as.numeric(max_vals[1,]), type = 'o')
  # lines(as.numeric(gbm_vals[1,]), type = 'o', col = 'red')

  # # density plots
  # plot(density(as.numeric(max_vals[1,])))
  # plot(density(as.numeric(gbm_vals[2,])), col = 'red')

  # loop through each model run
  for(m in 1:3){  # set m = 1 for debugging 
    # select model type
    m_name<-m_type[m]
    
    # select species by model type
    m_sp<-all_sp[[m]]
    
    # set up data variables 
    max_vals<-m_sp[,max_names]
    gbm_vals<-m_sp[,gbm_names]
    
    # # check number of bins
    # sqrt(length(max_vals))
    # sqrt(length(gbm_vals))
    
    # create image file
    png(paste0(outDir, "all_VariImp_curves_by_mod/", sp_nm[r], "_", m_name, "_all_VariImp.png"),
        width = 880, height = 880, units = "px", pointsize = 16)
    # set plot
    par(mfrow = c(2, 2))
    # loop through for each variable
    for(v in 1:4){  # set v = 1 for debugging
      # select variable
      v_name<-var_nm[v]
      
      # select by variable
      x_max<-as.numeric(max_vals[v,])
      max_mean<-mean(x_max)
      x_gbm<-as.numeric(gbm_vals[v,])
      gbm_mean<-mean(x_gbm)
      
      h_max<-hist(x_max, breaks = 10, plot = F)
      xfit_max<-seq(min(x_max), max(x_max), length = 40)
      yfit_max<-dnorm(xfit_max, mean = mean(x_max), sd = sd(x_max))
      yfit_max<-yfit_max*diff(h_max$mids[1:2])*length(x_max)
      # lines(xfit_max, yfit_max, col = "blue", lwd = 2)
      
      h_gbm<-hist(x_gbm, breaks = 10, plot = F)  
      xfit_gbm<-seq(min(x_gbm), max(x_gbm), length = 40)
      yfit_gbm<-dnorm(xfit_gbm, mean = mean(x_gbm), sd = sd(x_gbm))
      yfit_gbm<-yfit_gbm*diff(h_gbm$mids[1:2])*length(x_gbm)
      # lines(xfit_gbm, yfit_gbm, col = "blue", lwd = 2)
      
      # set lims
      d_xlim<-c(min(xfit_max, xfit_gbm), max(xfit_max, xfit_gbm))
      d_ylim<-c(min(yfit_max, yfit_gbm), max(yfit_max, yfit_gbm))
      
      # normal curve plots
      plot(xfit_max, yfit_max, xlim = d_xlim, ylim = d_ylim,
           col = "blue", lwd = max_mean*10, type = 'l', #lwd = 1.5, type = 'o',
           xlab = "Distribution", ylab = "Frequency",
           main = paste(sub("_", " ", sp_nm[r]), m_name, v_name, "Importance"))
      lines(xfit_gbm, yfit_gbm, col = "red", lwd = gbm_mean*10, type = 'l') #lwd = 1.5, type = 'o',)
      # legend("topleft", legend = c("MaxEnt", "GBM"), 
      #        col = c("blue", "red"), pch = 1, bty = "n")
      legend("bottomright", legend = c("MaxEnt", "GBM"), horiz = TRUE, 
             fill = c("blue", "red"), bty = "n", inset = c(0, 0.99), xpd = TRUE)
      legend("top", text.font = 4, bty = "n", inset = c(0, 1.08), xpd = T, 
             paste("MaxEnt Mean =", round(max_mean, 3), "GBM Mean =", round(gbm_mean, 3)))
    }
    # save plot
    dev.off()
  }
  rm(m, v)
  
  ### SAME ABOVE BUT BY BIO INSTEAD OF PER MODEL ###

  # set model colors
  m_col<-c("darkgreen", "red", "blue")
  
  for(v in 1:4){  # set v = 1 for debugging
    # select variable
    v_name<-var_nm[v]
    
    # create image file
    png(paste0(outDir, "all_VariImp_curves_by_bio/", sp_nm[r], "_", v_name, "_all_VariImp.png"),
        width = 680, height = 480, units = "px", pointsize = 12)
    # set plot
    par(mfrow = c(1, 2))
    
    # maxent
    for(m in 1:3){  # set m = 1 for debugging 
      # select model type
      m_name<-m_type[m]
      
      # select by variable
      x_max<-as.numeric(all_sp[[m]][v, max_names])
      max_mean<-mean(x_max)
      
      h_max<-hist(x_max, breaks = 10, plot = F)
      xfit_max<-seq(min(x_max), max(x_max), length = 40)
      yfit_max<-dnorm(xfit_max, mean = mean(x_max), sd = sd(x_max))
      yfit_max<-yfit_max*diff(h_max$mids[1:2])*length(x_max)
      # lines(xfit_max, yfit_max, col = "blue", lwd = 2)
      
      # # set lims for all models
      # d_xlim<-c(min(xfit_max, xfit_gbm), max(xfit_max, xfit_gbm))
      # d_ylim<-c(min(yfit_max, yfit_gbm), max(yfit_max, yfit_gbm))
      
      if(m == 1){
        par(mfrow = c(1, 2))
        # mean table
        all_max_mean<-max_mean
        # normal curve plots
        plot(xfit_max, yfit_max, xlim = c(0, 1), ylim = c(0, 45),
             col = m_col[m], lwd = max_mean*10, #lwd = 1.5, 
             type = 'l', xlab = "MaxEnt Models", ylab = "Frequency",
             main = paste(sub("_", " ", sp_nm[r]), v_name, "Importance"))
        legend("bottomright", legend = m_type, horiz = TRUE, 
               fill = c("darkgreen", "red", "blue"), bty = "n",
               inset = c(0, 0.99), xpd = TRUE)
      }else{
        # mean table
        all_max_mean<-c(all_max_mean, max_mean)
        # add lines to plot
        lines(xfit_max, yfit_max, col = m_col[m], 
              lwd = max_mean*10, type = 'l') #lwd = 1.5, type = 'o')
      }
    } # end maxent models
    # add mean variable importance
    legend("topright", paste(m_type, "mean =", round(all_max_mean, 3)), bty = "n")
    
    # gbm
    for(m in 1:3){  # set m = 1 for debugging 
      # select model type
      m_name<-m_type[m]
      
      # select by variable
      x_gbm<-as.numeric(all_sp[[m]][v,gbm_names])
      gbm_mean<-mean(x_gbm)
      
      h_gbm<-hist(x_gbm, breaks = 10, plot = F)  
      xfit_gbm<-seq(min(x_gbm), max(x_gbm), length = 40)
      yfit_gbm<-dnorm(xfit_gbm, mean = mean(x_gbm), sd = sd(x_gbm))
      yfit_gbm<-yfit_gbm*diff(h_gbm$mids[1:2])*length(x_gbm)
      # lines(xfit_gbm, yfit_gbm, col = "blue", lwd = 2)
      
      # # set lims for all models
      # d_xlim<-c(min(xfit_max, xfit_gbm), max(xfit_max, xfit_gbm))
      # d_ylim<-c(min(yfit_max, yfit_gbm), max(yfit_max, yfit_gbm))
      
      if(m == 1){
        # mean table
        all_gbm_mean<-gbm_mean
        # normal curve plots
        plot(xfit_gbm, yfit_gbm, xlim = c(0, 1), ylim = c(0, 45),
             col = m_col[m], lwd = gbm_mean*10, #lwd = 1.5, 
             type = 'l', xlab = "GBM Models", ylab = "Frequency",
             main = paste(sub("_", " ", sp_nm[r]), v_name, "Importance"))
        legend("bottomright", legend = m_type, horiz = TRUE, 
               fill = c("darkgreen", "red", "blue"), bty = "n",
               inset = c(0, 0.99), xpd = TRUE)
      }else{
        # mean table
        all_gbm_mean<-c(all_gbm_mean, gbm_mean)
        # add lines to plot
        lines(xfit_gbm, yfit_gbm, col = m_col[m], 
              lwd = max_mean*10, type = 'l') #lwd = 1.5, type = 'o')
      }
      
    } # end gbm models
    # add mean variable importance
    legend("topright", paste(m_type, "mean =", round(all_gbm_mean, 3)), bty = "n")
    
    # save plot
    dev.off()
  }
  
}
