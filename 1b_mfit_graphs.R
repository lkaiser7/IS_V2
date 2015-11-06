### graphs from model fitting ###
### use after 1_model_fitting ### 

##################
##### SET UP #####
##################

# load necessary packages
library(ggplot2)
library(reshape2)
library(plyr)

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
  evalMat0 = read.csv(paste0(project_path, 'all_eval_mat_', eval_stat, '.csv'))
  # loop through all models 
  for (model in models_to_run){
    # select evaluation statistic values per model
    evalMat = evalMat0[evalMat0[,2] == model,3:ncol(evalMat0)]
    # transpose data data frame
    evalMat = as.data.frame(t(evalMat))
    # rename for species runs
    names(evalMat) = evalMat0[evalMat0[,2]==model,1]
    # create two columns of species name and statistic values
    long = melt(evalMat)
    # rename new columns
    names(long) = c("Species", "Value")
    # find species place in vector of species names
    long = long[long[,1] %in% all_sp_nm,]
    # reorder species lists by name
    long = long[order(long[,1]),]
    # sort long data using species as factor levels
    long$Species<-factor(long$Species, levels = sort(levels(long$Species)))
    
    # name output image file
    tiff_name = paste0(outDir ,eval_stat, "_", model, "_variable_importance_box_plot.tif")
    # create blank image file
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # store basic qplot boxplot
    a = qplot(Species, Value, data = long, geom = c("boxplot"), fill = Species, main = "",
              xlab = "", ylab = paste(eval_stat, "model evaluation")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    # remove legend to plot
    a = a + theme(legend.position = "none")
    # save image file output
    ggsave(filename = tiff_name, plot = a)    
    
    # sign-posting of completed box plot for each evaluation statistic and model
    cat('done with', eval_stat, model, 'variable importance box plot \n')
  }
}

# load and format data in loop as done previously done above for violin plots
for (eval_stat in eval_stats){
  evalMat0=read.csv(paste0('all_eval_mat_',eval_stat,'.csv'))
  for (model in models_to_run){
    evalMat=evalMat0[evalMat0[,2]==model,3:ncol(evalMat0)]
    evalMat=as.data.frame(t(evalMat))
    names(evalMat)=evalMat0[evalMat0[,2]==model,1]
    long=melt(evalMat)
    names(long)=c("Species", "Value")
    long=long[long[,1] %in% spp_nm,]
    long=long[order(long[,1]),]
    long$Species <- factor(long$Species, levels = sort(levels(long$Species)))
    
    # name output image file
    tiff_name = paste0(outDir ,eval_stat, "_", model, "_variable_importance_violin_plot.tif")
    # create blank image file
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # store violin qplot
    a = qplot(Species, Value, data = long, geom = c("violin"), fill = Species, main="",
            xlab = "", ylab = paste(eval_stat, "model evaluation"))
    # + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    # add scale
    a = a + geom_violin(scale = "width")
    # add legend
    a = a + guides(fill = guide_legend(keywidth = 1, keyheight = 1.2))
    # remove axis elelments
    a = a + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    # rotate graph horizontally
    a = a + coord_flip()
    # resize plot and add corresponding legend
    a = a + guides(fill = guide_legend(reverse = TRUE), guide = guide_legend(title = NULL))
    # plot final graph with all previously stored adjustments
    a
    # save image file output
    ggsave(filename = tiff_name, plot = a)    
    
    # sign-posting of completed violin plot for each evaluation statistic and model
    cat('done with', eval_stat, model, 'variable importance box plot \n')
  }
}

###############################
##### VARIABLE IMPORTANCE #####
###############################

# open variable importance file created in script 1a
varImp0 = read.csv(paste0(project_path, 'all_VariImp.csv'))

# get column names of variable importance data
colnames = names(varImp0)
# create an empty list 
model_cols = list()
# loop through models run and add to list
for (model in models_to_run){ # set model = "GBM" for debugging
  # create temporary vector of number of column per model
  jnk = grep(paste0(model, "+"), colnames, perl = TRUE, value = FALSE)  
  # add model to list
  model_cols[[length(model_cols) + 1]] = jnk
}

# set species counter to 1
sp_nm = all_sp_nm[1]

# start looping though all species
for (sp_nm in all_sp_nm){
  # loop through all models run
  for (model in models_to_run){
    # number of models to loop through
    model_n = which(models_to_run==model)
    # select variable importance data from specific model
    varImp = varImp0[varImp0[,1]==sp_nm,model_cols[[model_n]]]
    # transform the data into a data frame
    varImp = as.data.frame(t(varImp))
    # rename columns from names stored above per bioclimatic variable
    names(varImp) = varImp0[varImp0[,1]==sp_nm, "rownames.Spp_VariImp."]
    # create two columns of bioclimatic predictors and their values
    long = melt(varImp)
    # rename new columns
    names(long)=c("Predictor", "Value")
    
    # name output image file
    tiff_name=paste0(outDir, sp_nm, "_",model, "_variable_importance_box_plot.tif")
    # create blank image file
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # store basic qplot boxplot
    a = qplot(Predictor, Value, data = long, geom = c("boxplot", "jitter"), 
            fill = Predictor, main = "", xlab = "", ylab = "Variable importance" )
    # save image file output
    ggsave(filename = tiff_name, plot = a)    
    
    # sign-posting of completed box plot for bioclimatic variables    
    cat('done with', sp_nm, model, 'variable importance box plot \n')
  }
}

# load and format data in loop as done previously done above for violin plot
for (sp_nm in spp_nm){
  for (model in models_to_run){
    model_n = which(models_to_run==model)
    varImp = varImp0[varImp0[,1]==sp_nm,model_cols[[model_n]]]
    varImp = as.data.frame(t(varImp))
    names(varImp) = varImp0[varImp0[,1]==sp_nm,"rownames.Spp_VariImp."]
    long = melt(varImp)
    names(long) = c("Predictor", "Value")
    
    # name output image file
    tiff_name=paste0(outDir, sp_nm, "_", model, "_variable_importance_violin_plot.tif")
    # create blank image file
    tiff(tiff_name, res = 300, units = "in", pointsize = 12,
         width = 10, height = 10, compression = "lzw")
    # store violin qplot
    a = qplot(Predictor, Value, data = long, geom = c("violin"),
            fill = Predictor, main = "", xlab = "", ylab = "Variable importance" )
    # add scale
    a = a + geom_violin(scale = "width")
    # remove axis elements 
    a = a+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    # save image file output
    ggsave(filename = tiff_name, plot = a) 

    # sign-posting of completed violin plot for bioclimatic variables    
    cat('done with', sp_nm, model, 'variable importance violin plot \n')
  }
}

####################################
##### END MODEL FITTING GRAPHS #####
####################################