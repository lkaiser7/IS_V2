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
  # replace species naming convention of "_" with "." 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
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
  
  # if only running one model, use this code for evaluation statistics
  if (length(models_to_run) == 1){
    # loop to combine all evaluations statistics from each species 
    for (eval_stat in eval_stats){ # set eval_stat = eval_stats for debugging
      # save all eval_stat values from model evaluation as a single vector
      Spp_eval<-c(myBiomodModelEval[eval_stat, "Testing.data",,,])
      Spp_eval<-rbind(Spp_eval)
      # combine species evaluation statistics with species name and statistic
      Spp_eval=cbind(sp_nm, models_to_run, Spp_eval)
      
      # get row (RUN) and column (PA) names of model evaluations
      pas<-colnames(myBiomodModelEval[eval_stat, "Testing.data",,,])
      runs<-rownames(myBiomodModelEval[eval_stat, "Testing.data",,,])
      # create column names for model evaluations
      for (runPA in 1:length(pas)){
        mat_cols<-paste(runs, pas[runPA], sep = ".")
        if (runPA == 1){
          matrix_cols<-mat_cols
        }else{
          matrix_cols<-c(matrix_cols, mat_cols)  
        }
      }
      # add species and run names to species evalutations
      colnames(Spp_eval)<-c("sp_nm", "models", matrix_cols)
      
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
      # store formatted evalutation statistics as a data frame
      
    }
  }else{
    # otherwise loop to combine all evaluations statistics from each model per species 
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
  }  
  
  ###############################
  ##### VARIABLE IMPORTANCE #####
  ###############################
  
  # run for each model 
  
  # save all variable importance values from models as a data frame
  Spp_VariImp<-data.frame(get_variables_importance(myBiomodModelOut))
  # sort columns in alphabetical order to keep models together
  Spp_VariImp<-Spp_VariImp[,order(names(Spp_VariImp))]
  
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

#####################
##### GET MEANS #####
#####################

# loop through and save each evaluation statistic type
for (eval_stat in eval_stats){
  # create file name for all model evaluations
  FileName<-paste0(project_path, "all_eval_mat_", eval_stat, ".csv")
  # save model evaluations per statistic as csv file
  write.csv(get(paste0("all_eval_mat_", eval_stat)), file = FileName, row.names = FALSE)
  
  # create temporary file for evalutation statistic
  tmp_eval_map = get(paste0("all_eval_mat_", eval_stat))
  
  # run only if models_to_run = 1 and matrix was created 
  if (class(tmp_eval_map) == "matrix"){
    # convert first two columns to a data frame
    tmp_df<-data.frame(tmp_eval_map[, 1:2])
    # convert matrix to numeric class
    class(tmp_eval_map)<-"numeric"
    # create data frame of all evaluation statistics 
    tmp_eval_map<-cbind(tmp_df, tmp_eval_map[, 3:dim(tmp_eval_map)[2]])
  }
  
  # rename first two columns
  names(tmp_eval_map)[c(1:2)] = c("species","model")
  # select first column for species names
  tmp_eval_sp = tmp_eval_map[, 1:2]
  
  # calculate mean per species by each model 
  model_means<-rowMeans(tmp_eval_map[, 3:dim(tmp_eval_map)[2]])
  # add mean to created data frame per species
  tmp_eval_byModel<-data.frame(species = tmp_eval_sp[,1], 
                               model = tmp_eval_sp[,2], model_means)
  
  # create file name for mean evaluations of statistics
  FileName<-paste0(project_path, "all_eval_meanBYmodel_mat_", eval_stat, ".csv")
  # save mean evaluation statistics as csv file
  write.csv(tmp_eval_byModel, file = FileName, row.names = FALSE)
  
  # reshape data frame of evaluation statistics per species
  tmp_eval_map<-reshape(tmp_eval_map, timevar=c("model"), idvar=c("species"), dir="wide")
  # get mean of evaluation statistic 
  meanEval = apply(tmp_eval_map[, 2:dim(tmp_eval_map)[2]], 1, mean, na.rm = TRUE)
  # add mean to created data frame per species
  tmp_eval_map2 = data.frame(species = tmp_eval_sp[,1], meanEval)
  # remove row assigned names
  row.names(tmp_eval_map2)<-NULL 
  
  # create file name for mean evaluations of statistics
  FileName<-paste0(project_path, "all_eval_mean_mat_", eval_stat, ".csv")
  # save mean evaluation statistics as csv file
  write.csv(tmp_eval_map2, file = FileName, row.names = FALSE)
}

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

####################################
##### END MODEL FITTING TABLES #####
####################################