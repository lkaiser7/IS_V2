### BIOMOD2: calculate shift ###
### see source script inputs ###
### step 6:

##################
##### SET UP #####
##################

# load necessary packages 
library(biomod2)
library(stringr)

# create directory for tables in project folder 
dir.create(paste0(outDir, "tables/"), showWarnings = F)

# set climate data directory to baseline bioclims
clim_data_dir = biobaseRun 

# set model resolution pixel area
px_area = model_resolution^2 
# store vector of variable to include in calculations
vars = c("Eval stat", "Species", "baseline_area", "future_area", "% change", 
         "lost", "kept", "gained", 
         "baseline_suitability", "future_suitability", "% suitability change")
# create empty data frame to compiled modeling results
all_stats<-as.data.frame(setNames(replicate(length(vars), numeric(0), simplify = F), vars)) 

####################################
##### CALCULATE SPECIES RANGES #####
####################################

# set species name to first species
sp_nm = all_sp_nm[1]
# set ensemble evaluation statistic to first stat
eval_stat = spp_ensemble_eval_stats[1]

# loop through all ensemble evaluation statistics 
for (eval_stat in spp_ensemble_eval_stats){
  # loop through all species 
  for (sp_nm in all_sp_nm){
    # convert species name to character
    sp_nm = as.character(sp_nm)  
    # store species name for file names
    sp_nm0 = sp_nm
    # replace species naming convention of "_" with "." 
    sp_nm = str_replace_all(sp_nm,"_", ".")
    
    # print sign posting of ongoing range shift calculations
    cat('\n Calculating range shifts for', sp_nm0, '...')
    
    # generate random deviates per observations 
    a = round(1000*runif(1))
    
    # store file output name for raster
    response_raster_nm = paste0('output_rasters/main/', sp_nm0, "_response_zones_",
                                eval_stat, "_", spp_ensemble_type, "_", comp_projects[2])
    # set raster path name tiff file in project path 
    response_raster_nm = paste0(project_path, response_raster_nm,".tif")
    # load tiff file as raster
    response_raster = raster(response_raster_nm)
        
    # create temporary raster of all species range area
    jnk = response_raster > 0
    # compute zonal statistics 
    jnk = zonal(jnk, response_raster, fun = sum, na.rm = TRUE)
    # store zonal statistics as a matrix
    jnk = matrix(jnk[jnk[, 1] > 0,], ncol = 2)
    # multiple zonal statistics value by pixel area
    jnk[,2] = jnk[, 2]*px_area
    # create empty data frame for zone area
    zone_area = as.data.frame(matrix(c(1, 2, 3, 0, 0, 0), nrow = 3, ncol = 2))
    # name column headers
    names(zone_area)=c("zone", "area")
    # store calculated area of zonal statistics values
    zone_area[jnk[, 1], 2] = jnk[, 2]
    # select column of only zonal area
    zone_area = zone_area[, 2]
    # calculate current zonal area
    curr = zone_area[1] + zone_area[2]
    # calculate future zonal area
    fut = zone_area[2] + zone_area[3]
    # combine calculated results (current, future, change, and zone area)
    results = c(curr, fut, (fut-curr)*100/curr, zone_area)
    
    # load raster mask of current range suitability
    baseline_masked_suitability = raster(paste0(project_path, 'output_rasters/', sp_nm0,"_",
                                                "clipped_suitability_",comp_projects[1],"_",
                                                eval_stat,"_",spp_ensemble_type, ".tif"))
    # load raster mask of future range suitability
    future_masked_suitability = raster(paste0(project_path, 'output_rasters/', sp_nm0,"_", 
                                              "clipped_suitability_",comp_projects[2],"_",
                                              eval_stat,"_",spp_ensemble_type, ".tif"))
    
    # create mask of current and future suitability
    baseline_masked_suitability[baseline_masked_suitability == 0] = NA
    future_masked_suitability[future_masked_suitability == 0] = NA    
    # calculate mean suitability score from raster masks 
    jnk1 = cellStats(baseline_masked_suitability, stat='mean', na.rm = TRUE)
    jnk2 = cellStats(future_masked_suitability, stat='mean', na.rm = TRUE)
    
    # store results in a matrix
    results = matrix(c(results, jnk1, jnk2, 100*(jnk2-jnk1)/jnk1), nrow = 1)
    # convert matrix to data frame
    results = as.data.frame(results)
    # add evaluation statistic and species name
    results = cbind(eval_stat, sp_nm0, results)
    # name column headers 
    names(results)=vars
    
    # store file output name for csv table
    out_csv_name = paste0(project_path, 'tables/', sp_nm0, "_metrics_", eval_stat, "_", 
                          spp_ensemble_type, "_", comp_projects[2], ".csv")
    # save table to file
    write.table(results, file = out_csv_name, sep=",", row.names = FALSE)
    
    # add results to file with all species and evaluation statistic calculations
    all_stats = rbind(all_stats, results)
  }
}

# store file output name for all calculated statistics per species and metrics
out_csv_name = paste0(project_path, 'tables/', "all_spp_distr_shift_metrics_", 
                      comp_projects[2], ".csv")
# save table to file
write.table(all_stats, file = out_csv_name, sep = ",", row.names = FALSE)

########################################
##### END RANGE SHIFT CALCULATIONS #####
########################################