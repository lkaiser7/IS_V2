### BIOMOD2: creates rasters ###
### see source script inputs ###
### step 4:

##################
##### SET UP #####
##################

# load necessary packages 
library(biomod2)
library(stringr)
library(colorRamps)
library(raster)
######################
##### FUNCTIONS ######
######################

#########################
##### RUN FUNCTIONS #####
#########################


# set species name to first species
sp_nm = all_sp_nm[1]

# select first evaluation statistic
eval_stat = spp_ensemble_eval_stats[1]
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
    cat('\n',sp_nm,'grid conversion')
    
    
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
      
      #file_name=file_name1
      replace_grd_fx=function(file_name){
        if (file.exists(file_name)){
          library(terra) #using terra package to save as raster package does not keep layer names
          #https://stackoverflow.com/questions/26763013/r-write-rasterstack-and-preserve-layer-names
          # store first file name for raster .grd file
          temp_raster = rast(file_name) 
          output_nm=sub(pattern = ".grd$", replacement = ".tif", file_name)
          if (file.exists(output_nm) == FALSE | overwrite == 1){
            terra::writeRaster(temp_raster, filename = output_nm, gdal=c(compress="LZW"), overwrite=T)
          }
          unlink(file_name, recursive=T, force=T)
          unlink(sub(pattern = ".grd$", replacement = ".gri", file_name), recursive=T, force=T)
        }
      }
      
      file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                          proj_nm, "_", sp_nm, "_ensemble.grd")
      replace_grd_fx(file_name1)
      
      file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                          proj_nm, "_", sp_nm, "_ensemble_", eval_stat, "bin.grd")
      replace_grd_fx(file_name1)
      
      ####delete unused grids
      file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                          proj_nm, "_", sp_nm, "_", eval_stat, "bin.grd")
      unlink(file_name1, recursive=T, force=T)
      unlink(sub(pattern = ".grd$", replacement = ".gri", file_name1), recursive=T, force=T)
      
      file_name1 = paste0(project_run, "/", sp_nm, "/proj_", proj_nm, "/proj_", 
                          proj_nm, "_", sp_nm, ".grd")
      unlink(file_name1, recursive=T, force=T)
      unlink(sub(pattern = ".grd$", replacement = ".gri", file_name1), recursive=T, force=T)
      
    }
  }
}

##############################
##### END RASTER OUTPUTS #####
##############################