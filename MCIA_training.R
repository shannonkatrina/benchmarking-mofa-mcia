#############################
##   MCIA model training   ##
##     R version 4.0.2     ##
#############################

###############
## libraries ##
###############

library(tidyverse)
library(omicade4)
library(rlist)
library(ggplot2)
library(ggpubr)
source("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/MODI_helper_functions.R")

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

# convert to list of dfs
MCIA_arrays <- tidy_mae(mae_subset)
sapply(MCIA_arrays, dim) # check dimensions of list

########################
##  Training a model  ##
########################

# run MCIA
mcoin <- mcia(MCIA_arrays, cia.nf=10) # include 10 factors for analysis

# save
write_rds(mcoin,"/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/MCOIN.rds")


