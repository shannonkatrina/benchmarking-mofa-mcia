#############################
##   MOFA model training   ##
##     R version 4.0.2     ##
#############################

###############
## libraries ##
###############

library(ggplot2)
library(tidyverse)
library(MOFA2)
library(rlist)
library(reticulate)
source("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/MODI_helper_functions.R")
use_condaenv("/Users/shannonbrown/opt/anaconda3/envs/mofa2/", required=TRUE) 
py_config()

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

###############
## load data ##
###############

mae <- readRDS("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/mae_all_2022.08.05_naomit.rds")
mae@ExperimentList

# filter top variable features in mae 
mae <- mae[,,c("ATAC","EP","PP","RNA")] # remove metabolomics
top_varFeats <- list("ATAC" = top_varFeats_view(mae, "ATAC", NA_omit=TRUE, 5000),
                     "EP" = top_varFeats_view(mae, "EP", NA_omit=TRUE, 2000),
                     #"Met" = top_varFeats_view(mae, "Met", NA_omit=TRUE, 1000),
                     "PP" = top_varFeats_view(mae, "PP", NA_omit=TRUE, 2000),
                     "RNA" = top_varFeats_view(mae, "RNA", NA_omit=TRUE, 5000))

# get list of cell lines in all platforms
intersect_cells <- Reduce(intersect, lapply(names(mae), function(x){colnames(assay(mae[[x]]))})) 

# filter multiassay experiment
mae_subset <- mae[top_varFeats, intersect_cells]
mae_subset@ExperimentList

########################
##  Training a model  ##
########################

# create MOFA object
MOFAobject <- create_mofa_from_MultiAssayExperiment(mae_subset)
views_names(MOFAobject)

plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 11
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$seed <- 1
train_opts$drop_factor_threshold <- 0 # removes factor if explains less than x variance
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

MOFAobject.trained <- run_mofa(MOFAobject)

# save
saveRDS(MOFAobject.trained, "/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/MOFAtrained.rds")


