### graphs from model fitting ###
### use after 1_model_fitting ### 

##################
##### SET UP #####
##################

# load necessary packages
library(ggplot2)
library(reshape2)
library(plyr)

# create folder in output directory for response curve plots
vi_fold<-paste0(outDir, "VarImp_plots/")
dir.create(vi_fold, showWarnings = FALSE) 

#################################
##### EVALUATION STATISTICS #####
#################################

# select first evaluation statistic
eval_stat = eval_stats[1]
# select first model from runs
model = models_to_run[1]

# loop through all evaluation statistics for box plots
for (eval_stat in eval_stats){
  # open evaluation statistic file created in script 1a
  evalMat0 = read.csv(paste0(outDir, 'all_eval_mat_', eval_stat, '.csv'))
  # loop through all models 
  for (model in models_to_run){
    # select evaluation statistic values per model
    evalMat = evalMat0[evalMat0[,3] == model,4:ncol(evalMat0)]
    # transpose data data frame
    evalMat = as.data.frame(t(evalMat))
    # rename for species runs
    names(evalMat) = evalMat0[evalMat0[,3] == model,1]
    # create two columns of species name and statistic values
    long = melt(evalMat)
    # rename new columns
    names(long) = c("Species", "Value")

    # name output image file (ALL GGPLOT FILES MUST BE SAVED AS '.tiff' FILES!)
    t_name = paste0(vi_fold, eval_stat, "_", model, "_model_skil_box_plot.tiff")
    # create blank image file
    # tiff(t_name, res = 300, units = "in", pointsize = 12,
    #      width = 10, height = 10, compression = "lzw")
    # store basic qplot boxplot
    a = qplot(Species, Value, data = long, geom = c("boxplot"), fill = Species, main = "",
              xlab = "", ylab = paste(eval_stat, "model evaluation")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    # remove legend to plot
    a = a + theme(legend.position = "none")
    #print(a)
    
    # save image file output
    ggsave(filename = t_name, plot = a)    
    graphics.off()
    
    # sign-posting of completed box plot for each evaluation statistic and model
    cat('done with', eval_stat, model, 'model_skil box plot \n')
    
    if (model == models_to_run[1]){
      all_model_long=long
    }else{
      all_model_long=rbind(all_model_long, long)
    }
  }
  # name output image file (ALL GGPLOT FILES MUST BE SAVED AS '.tiff' FILES!)
  t_name = paste0(vi_fold, eval_stat, "_allModels_model_skil_box_plot.tiff")
  # create blank image file
  # tiff(t_name, res = 300, units = "in", pointsize = 12,
  #      width = 10, height = 10, compression = "lzw")
  # store basic qplot boxplot
  a = qplot(Species, Value, data = all_model_long, geom = c("boxplot"), fill = Species, main = "",
            xlab = "", ylab = paste(eval_stat, "model evaluation")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  # remove legend to plot
  a = a + theme(legend.position = "none")
  a
  #print(a)
  
  # save image file output
  ggsave(filename = t_name, plot = a)    
  graphics.off()
  
}


###############################
##### VARIABLE IMPORTANCE #####
###############################

# open variable importance file created in script 1a
varImp0 = read.csv(paste0(outDir, 'all_VariImp.csv'))

# get column names of variable importance data
colnames = names(varImp0)
# create an empty list 
model_cols = list()
# loop through models run and add to list
model = models_to_run[1]
for (model in models_to_run){ # set model = "GBM" for debugging
  # create temporary vector of number of column per model
  jnk = grep(paste0(model, "+"), colnames, perl = TRUE, value = FALSE)  
  # add model to list
  model_cols[[length(model_cols) + 1]] = jnk
}

# set species counter to 1
sp_nm = all_sp_nm[1]
model = models_to_run[1]

# start looping though all species
for (sp_nm in all_sp_nm){
  # loop through all models run
  for (model in models_to_run){
    # number of models to loop through
    model_n = which(models_to_run==model)
    # select variable importance data from specific model
    varImp = varImp0[varImp0[,2]==sp_nm,model_cols[[model_n]]]
    # transform the data into a data frame
    varImp = as.data.frame(t(varImp))
    # rename columns from names stored above per bioclimatic variable
    names(varImp) = varImp0[varImp0[,2]==sp_nm, "rownames.Spp_VariImp."]
    # create two columns of bioclimatic predictors and their values
    long = melt(varImp)
    # rename new columns
    names(long)=c("Predictor", "Value")
    
    # name output image file
    t_name=paste0(vi_fold, sp_nm, "_",model, "_variable_importance_box_plot.tiff")
    a = qplot(Predictor, Value, data = long, geom = c("boxplot", "jitter"), 
            fill = Predictor, main = "", xlab = "", ylab = "Variable importance" )
    ggsave(filename = t_name, plot = a, dpi = 300, units = "in", width = 10, height = 10, compression = "lzw")    
    graphics.off()
    
    # sign-posting of completed box plot for bioclimatic variables    
    cat('done with', sp_nm, model, 'variable importance box plot \n')
  }
}


####################################
##### END MODEL FITTING GRAPHS #####
####################################