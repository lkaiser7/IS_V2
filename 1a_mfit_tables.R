### tables from model fitting ###
### use after 1_model_fitting ### 

##################
##### SET UP #####
##################

# load necessary packages
library(biomod2)
library(reshape2)
library(plyr)

# select individula species (miconia = all_sp_nm[7])
sp_nm = all_sp_nm[1]
# set counter (i) to 1 to start
i = 1
# begin loop through all species names
for (sp_nm in all_sp_nm){
  # set name of file saved from model fitting script 
  workspace_name = paste0(project_path, sp_dir, sp_nm, "_modelfitting.RData") 
  # load workspace data from model runs
  load(workspace_name)
  # run evaluations on BIOMOD2 model output
  myBiomodModelEval<-get_evaluations(myBiomodModelOut)   
  # list dimensions of model evaluations
  dimnames(myBiomodModelEval)
  
  #################################
  ##### EVALUATION STATISTICS #####
  #################################
  
  # loop to combine all evaluations statistics from each species 
  for (eval_stat in eval_stats){ # set eval_stat = eval_stats for debugging
    # save all eval_stat values from model evaluation as a data frame
    Spp_eval<-data.frame(myBiomodModelEval[eval_stat, "Testing.data",,,])
    # add species and run names to data frame
    Spp_eval=cbind(matrix(sp_nm, dim(Spp_eval)[1],1), rownames(Spp_eval), Spp_eval)
    
    # create master file for evaluations from the first species
    if (i==1){
      # give new name to data frame for master file
      assign(paste0("all_eval_mat_", eval_stat), Spp_eval)
      # otherwise for all other species
    }else{
      # create a temporary data frame with evaluation data
      jnk = rbind(get(paste0("all_eval_mat_", eval_stat)), Spp_eval)
      # add temporary data frame to master file of evaluations
      assign(paste0("all_eval_mat_", eval_stat), jnk)      
    }
  }
  
  ###############################
  ##### VARIABLE IMPORTANCE #####
  ###############################
  
  # save all variable importance values from models as a data frame
  Spp_VariImp<-data.frame(get_variables_importance(myBiomodModelOut))
  # add species and run names to data frame
  Spp_VariImp=cbind(matrix(sp_nm, dim(Spp_VariImp)[1],1), rownames(Spp_VariImp), Spp_VariImp)
  
  # create master file for variable importance from the first species
  if (i==1){
    # give new name to data frame for master file
    all_var_imp_mat = Spp_VariImp
    # otherwise for all other species
  }else{
    # add other variable importance values to master file
    all_var_imp_mat = rbind(all_var_imp_mat, Spp_VariImp)
  }
  
  # add 1 to counter (i) to go to next species in loop
  i = i + 1
}
# create file name for all model variable importance values
FileName<-paste0(project_path, "all_VariImp.csv")
# save variable importance as csv file
write.csv(all_var_imp_mat, file = FileName, row.names = FALSE)

# select columns of species names and bioclimatic variables
all_var_imp_mean = all_var_imp_mat[,1:2]
# get mean of variable importance values
meanVarImp = apply(all_var_imp_mat[, 3:dim(all_var_imp_mat)[2]],1, mean, na.rm = TRUE)
# add mean values to data frame of species names and bioclimatic variables
all_var_imp_mean = cbind(all_var_imp_mean, meanVarImp)
# rename the column headers
names(all_var_imp_mean) = c("species", "var", "meanVarImp")
# remove row assigned names
row.names(all_var_imp_mean)<-NULL 
# create file name for mean of variable importance values
FileName<-paste0(project_path, "all_VariImp_mean.csv")
# save variable importance means as csv file
write.csv(all_var_imp_mean, file = FileName, row.names = FALSE)

# loop through and save each evaluation statistic type
for (eval_stat in eval_stats){
  # create file name for all model evaluations
  FileName<-paste0(project_path, "all_eval_mat_", eval_stat, ".csv")
  # save model evaluations per statistic as csv file
  write.csv(get(paste0("all_eval_mat_", eval_stat)), file = FileName, row.names = FALSE)
  
  # create temporary file for evalutation statistic
  tmp_eval_map = get(paste0("all_eval_mat_", eval_stat))
  # rename first two columns
  names(tmp_eval_map)[c(1:2)] = c("species","model")
  # reshape data frame of evaluation statistics per species
  tmp_eval_map<-reshape(tmp_eval_map, timevar=c("model"), idvar=c("species"), dir="wide")
  # select first column for species names
  tmp_eval_map2 = tmp_eval_map[,1]
  # get mean of evaluation statistic 
  meanEval = apply(tmp_eval_map[, 2:dim(tmp_eval_map)[2]],1, mean, na.rm = TRUE)
  # add mean to created data frame per species
  tmp_eval_map2 = data.frame(species = tmp_eval_map2, meanEval)
  # remove row assigned names
  row.names(tmp_eval_map2)<-NULL 
  
  # create file name for mean evaluations of statistics
  FileName<-paste0(project_path, "all_eval_mean_mat_", eval_stat, ".csv")
  # save mean evaluation statistics as csv file
  write.csv(tmp_eval_map2, file = FileName, row.names = FALSE)
}

####################################
##### END MODEL FITTING TABLES #####
####################################