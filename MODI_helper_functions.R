# plot PCA
# matrux colnames must be stripped_cell_line_name and metadta must contain this field to match on

plot_PCA_from_matrix <- function(matrix, ntop, metadata, intgroup, title=NULL, subtitle=NULL, caption=NULL, shape_group=NA){
  
  # matrix must be in format rows = features, cols = samples
  # calculate the variance for each gene
  rv <- rowVars(matrix)
  
  # select the ntop features by variance
  top_feats <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data for the selected features
  pca <- prcomp(t(matrix[top_feats,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2 )
  
  # check metadata in same order as matrix
  metadata <- metadata[match(colnames(matrix), metadata$stripped_cell_line_name),]
  
  # pull intgroup
  group <- metadata %>%
    pull(intgroup)
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, name=colnames(matrix))
  
  # plot
  
  if(is.na(shape_group)){
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
      geom_point(size=3) + 
      xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
      labs(color=intgroup,
           title = title,
           subtitle = subtitle,
           caption = caption) 
  
  }else{
    
    shape = metadata %>%
      pull(shape_group)
      
    ggplot(data=d, aes_string(x="PC1", y="PC2", color="group", shape="shape")) + 
      geom_point(size=3) + 
      xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
      labs(color=intgroup,
           shape=shape_group,
           title = title,
           subtitle = subtitle,
           caption = caption) 
    
  }
  
}


# filter top variable features of multiassay experiment
# use when want to filter all views to the same number of features
# function creates list of top n variable features within each omics in the mae. This list can then be used to filter the mae
# NA_omit = logical option to filter out features that have NA values

find_top_varFeats <- function(mae, NA_omit, n_feats) {
  
  list <- vector(mode = "list")
  
  if(NA_omit){
    
    for (view in names(mae)){
      data <- assay(mae[[view]]) %>% 
        as.matrix() %>% 
        na.omit() %>% 
        matrixStats::rowVars(.) %>% 
        order(decreasing=T) %>% 
        head(n_feats) %>% 
        {rownames(mae[[view]])[.]} 
      
      list <- list.append(list,data)
    }
    
  }else{
    
    for (view in names(mae)){
      data <- assay(mae[[view]]) %>% 
        as.matrix() %>% 
        matrixStats::rowVars(.) %>% 
        order(decreasing=T) %>% 
        head(n_feats) %>% 
        {rownames(mae[[view]])[.]} 
      
      list <- list.append(list,data)
    }
    
  }
  
  names(list) <- names(mae)
  
  return(list)
  
}

#top_varFeats <- find_top_varFeats(mae, 1000)
#mae_subset <- mae[top_varFeats]

##########################################################################################

# filter top variable features of view within multiassay experiment
# useful when want to filter each view with different number of features
# function creates a vector of top n variable features. Combine into a list to then filter the mae
# NA_omit = logical option to filter out features that have NA values

top_varFeats_view <- function(mae, view, NA_omit, n_feats) {
  
  if(NA_omit){
    
    matrix <- na.omit(as.matrix(assay(mae[[view]])))
    
    topFeats_index <- matrix %>% 
      matrixStats::rowVars(.) %>%
      order(decreasing=T) %>%
      head(n_feats) 
      
    topFeats <- rownames(matrix[topFeats_index,])
    
  }else{
    
    matrix <- as.matrix(assay(mae[[view]]))
    
    topFeats_index <- matrix %>%
      matrixStats::rowVars(.) %>% 
      order(decreasing=T) %>% 
      head(n_feats) 
    
    topFeats <- rownames(matrix[topFeats_index,])
    
  }
  
  return(topFeats)
  
} 

#top_varFeats <- list("SE_ATACseq" = top_varFeats_view(mae, "SE_ATACseq", 5000),
#                     "SE_expressionProteomics" = top_varFeats_view(mae, "SE_expressionProteomics", 1000),
#                     "SE_metabolomics" = top_varFeats_view(mae, "SE_metabolomics", 1000),
#                     "SE_phosphoproteomics" = top_varFeats_view(mae, "SE_phosphoproteomics", 1000),
#                     "SE_transcriptomics" = top_varFeats_view(mae, "SE_transcriptomics", 5000))
#mae_subset <- mae[top_varFeats]

##########################################################################################

# function to remove features and cell lines with missing data
## removes features missing in certain cell lines
## removes cell lines missing omics layers
## outputs list of dataframes for MCIA

tidy_mae <- function(mae){
  
  list <- vector(mode = "list")
  
  # loop through each dataset in the mae and remove features with NA
  # ouputs list of data frames
  for (view in names(mae)){
    data <- assay(mae[[view]]) %>% as.matrix() %>% na.omit() %>% data.frame()
    list <- list.append(list,data)
  }
  
  names(list) <- names(mae) # name elements of list
  intersect_cells <- Reduce(intersect, lapply(list, colnames)) # get list of cell lines in all platforms
  
  list <- lapply(list, function(x) x[colnames(x) %in% intersect_cells]) # filter cell lines
  
  return(list)
}

##########################################################################################

# function to convert multiassay experiemnt to list of dataframes
mae_to_list <- function(mae){
  
  list <- vector(mode = "list")
  
  # loop through each dataset in the mae and remoVEe features with NA
  # ouputs list of data frames
  for (view in names(mae)){
    data <- assay(mae[[view]]) %>% data.frame()
    list <- list.append(list,data)
  }
  
  names(list) <- names(mae) # name elements of list
  
  return(list)
  
}  
