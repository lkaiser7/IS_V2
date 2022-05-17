project_dirs=c("C:/Users/lkaiser-local/Desktop/Phase1_SDMs/", "E:/invasives_SDM/")
rootDir=project_dirs[min(which(dir.exists(project_dirs)))]

# set working directory to main analysis folder
setwd(rootDir)

inat_df=read.csv("data/inat_data/inaturalist_gbif_simple.csv")
inat_sp_list=read.csv("data/inat_data/species_list.csv")

#View(inat_df)
inat_df=inat_df[!is.na(inat_df$xy_uncertainty_m),]
inat_df=inat_df[inat_df$xy_uncertainty_m<250,]
inat_df=inat_df[!is.na(inat_df$species),]
unique(inat_df$species)

to_replace=inat_sp_list$synonym
replacement=inat_sp_list$species
vector_to_apply=inat_df$species
#this works
#https://stackoverflow.com/questions/7547597/dictionary-style-replace-multiple-items
library(plyr)
foo <- mapvalues(vector_to_apply, from=to_replace, to=replacement)
unique(foo)

inat_df$species=foo

#subset by region

#create spatial dataframe
library(sf)
library(sp)
library(raster)
ac=inat_df[, -4]
names(ac)=c("species", "y", "x")
y=ac$y
x=ac$x
coordinates(ac) <- ~x+y
projection(ac) <- CRS('+proj=longlat +datum=WGS84')


#View(ac@data)

scopes=c("all_data", "hi_data", "no_hi_data")

hawaii_extent=extent(-160, -154, 18, 23)

crop(ac, hawaii_extent)

ac@data$species=sub(" ", "_", ac@data$species)

CP <- as(hawaii_extent, "SpatialPolygons")
proj4string(CP) <- CRS('+proj=longlat +datum=WGS84')

all_sp_nm = c('Clidemia_hirta', 'Falcataria_moluccana', 'Hedychium_gardnerianum',
              'Lantana_camara', 'Leucaena_leucocephala', 'Melinis_minutiflora',
              'Morella_faya', 'Panicum_maximum',
              'Passiflora_tarminiana', 'Pennisetum_clandestinum', 'Pennisetum_setaceum',
              'Psidium_cattleianum', 'Setaria_palmifolia','Schinus_terebinthifolius',
              'Cyathea_cooperi', 'Miconia_calvescens', 'Ulex_europaeus')

dir.create("data/inat_data/all_data", showWarnings = F)
dir.create("data/inat_data/hi_data", showWarnings = F)
dir.create("data/inat_data/no_hi_data", showWarnings = F)

sp_nm = all_sp_nm[16]
unique(ac@data$species)
for (sp_nm in all_sp_nm){
  cat("doing ", sp_nm, "\n")
  sp_ac=ac[ac@data$species==sp_nm,]
  
  if (nrow(sp_ac@data)>0){
    ## Clip the map
    library(rgeos)
    hi_data <- gIntersection(sp_ac, CP)#, byid=TRUE)
    #hi_data=as.spatialpo
    
    #gDifference from the rgeos
    no_hi_data <- gDifference(sp_ac, CP, byid=TRUE)
    
    ## Plot the output
    #plot(hi_data, col="khaki", bg="azure2")
    #plot(no_hi_data, col="khaki", bg="azure2")
    #plot(ac, col="khaki", bg="azure2")
    
    df=cbind(coordinates(hi_data), sp_nm)
    write.csv(df, paste0("data/inat_data/hi_data/",sp_nm,".csv"), row.names = F)
    
    df=cbind(coordinates(no_hi_data), sp_nm)
    write.csv(df, paste0("data/inat_data/no_hi_data/",sp_nm,".csv"), row.names = F)
    
    df=cbind(coordinates(sp_ac), sp_nm)
    write.csv(df, paste0("data/inat_data/all_data/",sp_nm,".csv"), row.names = F)
  }
}

# now merge
dir.create("data/merged_data/all_data", showWarnings = F, recursive = T)
dir.create("data/merged_data/hi_data", showWarnings = F, recursive = T)
dir.create("data/merged_data/no_hi_data", showWarnings = F, recursive = T)

sp_nm = all_sp_nm[2]
scope = scopes[3]


for (scope in scopes){
  for (sp_nm in all_sp_nm){
    cat("doing ", sp_nm, "\n")
    file_nm=paste0("data/orig_data/",scope,"/",sp_nm,".csv")
    if (file.exists(file_nm)){
      orig_df=read.csv(file_nm)
      orig_df=orig_df[,c(-5)]
      orig_exists=T
    }else{
      orig_exists=F
    }
    
    file_nm=paste0("data/inat_data/",scope,"/",sp_nm,".csv")
    if (file.exists(file_nm)) {
      inat_df=read.csv(file_nm)
      inat_df$temp=1
      inat_df=inat_df[,c(4, 3,1,2)]
      names(inat_df)=c("X", "Species", "decimalLongitude", "decimalLatitude")
      inat_exists=T
    }else{
      inat_exists=F
    }
    #View(orig_df)
    #View(inat_df)
    if (orig_exists&inat_exists){
      df=rbind(orig_df, inat_df)
    }else{
      if (orig_exists){
        df=orig_df
      }
      if (inat_exists){
        df=inat_df
      }
    }
    write.csv(df, paste0("data/merged_data/",scope,"/",sp_nm,".csv"), row.names = F)
  }
}
