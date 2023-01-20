###############################
## MAE creation from 5 omics ##
###############################

###########
## paths ##
###########

repo <- "/home/skb45807/MODI/modi_mres_project/"
data_dir <- "/home/skb45807/MODI/data/47_cancer_cells/"

###############
## libraries ##
###############

library(limma)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)
library(readr)
library(vsn)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(UpSetR)
source(paste0(repo,"MODI_helper_functions.R"))

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

###############
## load data ##
###############

metadata <- read.csv(paste0(data_dir,"new_cell_metadata_2022.08.03.csv"), na.strings=c("","NA"))

#####################
## transcriptomcis ##
#####################

# load data
log2TPM <- read.table(paste0(data_dir,"transcriptomics/nira_dxd_48celllines_log2tpmp1_baselineexpression_dmcbatchcorrected_filtered_updated.tsv"), header = T)
rownames(log2TPM) <- sub("([A-Z])0(.*)", "\\1\\2", gsub("_gsk", "", rownames(log2TPM))) # remove "gsk_" and any 0 in middle of string to improve metadata matching

# remove missing cell lines from df 
missing_celllines <- setdiff(rownames(log2TPM), metadata$stripped_cell_line_name)
log2TPM <- log2TPM[!rownames(log2TPM) %in% missing_celllines,]

# convert to matrix and remove NA
log2TPM_matrix <- t(na.omit(log2TPM))

#create se
metadata_transcriptomics <- column_to_rownames(filter(metadata, stripped_cell_line_name %in% colnames(log2TPM_matrix)),
                                               "stripped_cell_line_name")
se_transcriptomics <- SummarizedExperiment(log2TPM_matrix, colData = metadata_transcriptomics)

#################
## epigenomics ##
#################

# load data
ATACseq <- read.table(paste0(data_dir,"ATAC-seq/Combined_Filtered_BrCa_LuCa_OvCa_RawReads8_Deseq2LibSize_VST_Normalized_Median_Peak_Counts_LongFormat.tsv"), header = T, row.names = 1)

# remove missing cell lines from df 
missing_celllines <- setdiff(colnames(ATACseq), metadata$stripped_cell_line_name)
ATACseq <- ATACseq[,!colnames(ATACseq) %in% missing_celllines]

# filter promotor peaks 
promoter_peakIDs <- scan(paste0(data_dir,"ATAC-seq/promoter_PeakIDs.txt"), character()) # load id txt file
ATACseq <- ATACseq[rownames(ATACseq) %in% promoter_peakIDs, ]

# convert to matrix and remove NA
ATACseq_matrix <- as.matrix(na.omit(ATACseq))

# create se
metadata_epigenetics <- column_to_rownames(filter(metadata, stripped_cell_line_name %in% colnames(ATACseq_matrix)),
                                           "stripped_cell_line_name")
se_epigenetics <- SummarizedExperiment(ATACseq_matrix, colData = metadata_epigenetics)

###########################
## expression proetomics ##
###########################

# load data
load(paste0(data_dir,"expression_proteomics/EP_panel47_eprot_2021-06-25.RData"))
exp_proteomics <- EP47$data_long

# create average for cell line replicates
log2sia_ep <- exp_proteomics %>%
  select(stripped_cell_line_name_for_mapping, gene_name, log2_sia_bc) %>%
  filter(stripped_cell_line_name_for_mapping != "pooled") %>%
  group_by(stripped_cell_line_name_for_mapping, gene_name) %>%
  summarise(avg_sia = mean(log2_sia_bc)) %>%
  pivot_wider(names_from = stripped_cell_line_name_for_mapping, values_from = avg_sia) %>%
  column_to_rownames('gene_name')

# convert to matrix and remove NA
log2sia_ep_matrix <- as.matrix(na.omit(log2sia_ep))

# remove missing cell lines 
missing_celllines <- setdiff(colnames(log2sia_ep), metadata$stripped_cell_line_name)
log2sia_ep_matrix <- log2sia_ep_matrix[,!colnames(log2sia_ep_matrix) %in% missing_celllines]

# create se
metadata_EP <- column_to_rownames(filter(metadata, stripped_cell_line_name %in% colnames(log2sia_ep_matrix)),
                                  "stripped_cell_line_name")
se_EP <- SummarizedExperiment(log2sia_ep_matrix, colData = metadata_EP)

#######################
## phosphoproetomics ##
#######################

#  load data
pprot <- readRDS(paste0(data_dir,"phosphoproteomics/PP_panel47_pprot_2021-07-20_wo_outlier.rds"))
phosphoproteomics <- pprot$data_long

# create average for cell line replicates 
log2sia_pp <- phosphoproteomics %>%
  filter(!is.na(log2_sia_norm_aggr_bc)) %>%  # normalised & batch corrected
  filter(stripped_cell_line_name_for_mapping != "pooled") %>%
  mutate(psite_gene = paste(gene_name, psite, sep = '_')) %>%
  group_by(stripped_cell_line_name_for_mapping, psite_gene) %>%
  summarise(avg_sia = mean(log2_sia_norm_aggr_bc)) %>%
  pivot_wider(names_from = stripped_cell_line_name_for_mapping, values_from = avg_sia) %>%
  column_to_rownames('psite_gene')

# convert to matrix and remove NA
log2sia_pp_matrix <- as.matrix(na.omit(log2sia_pp))

# remove missing cell lines
missing_celllines <- setdiff(colnames(log2sia_pp), metadata$stripped_cell_line_name)
log2sia_pp_matrix <- log2sia_pp_matrix[,!colnames(log2sia_pp_matrix) %in% missing_celllines]

# create se
metadata_PP <- column_to_rownames(filter(metadata, stripped_cell_line_name %in% colnames(log2sia_pp_matrix)),
                                  "stripped_cell_line_name")
se_PP <- SummarizedExperiment(log2sia_pp_matrix, colData = metadata_PP)

##################
## metabolomics ##
##################

#  load data
log10intensity <- read.table(paste0(data_dir,"metabolomics/nira_dxd_metabolomics_cz_47celllines_baseline_log10intensity.tsv"), header = T, row.names = 1)
colnames(log10intensity) <- sub("([A-Z])0(.*)", "\\1\\2", colnames(log10intensity)) # fix "HCC0827"

# convert to matrix and remove NA
log10intensity_matrix <- as.matrix(na.omit(log10intensity))

# remove missing cell lines
missing_celllines <- setdiff(colnames(log10intensity), metadata$stripped_cell_line_name)
log10intensity_matrix <- log10intensity_matrix[,!colnames(log10intensity_matrix) %in% missing_celllines]

# create se
metadata_met <- column_to_rownames(filter(metadata, stripped_cell_line_name %in% colnames(log10intensity_matrix)),
                                   "stripped_cell_line_name")
se_met <- SummarizedExperiment(log10intensity_matrix, colData = metadata_met)

#########
## mae ##
#########

all_omics_se <- list(RNA=se_transcriptomics,
                     ATAC=se_epigenetics,
                     EP=se_EP,
                     PP=se_PP,
                     Met=se_met)

rownames(metadata) <- metadata$stripped_cell_line_name
mae <- MultiAssayExperiment(all_omics_se, colData = metadata)

colData(mae) # check rows and columns
mae@ExperimentList
saveRDS(mae, file = paste0(data_dir,"mae_all_2022.08.05_naomit.rds"))

################
## UpSet plot ##
################

joint_cell_ids <- lapply(all_omics_se, function(x){colnames(assay(x))})

jpeg(paste0(vis_dir,"se_intersection.jpeg"), width = 500, height=400)
upset(fromList(joint_cell_ids))
dev.off()

#################
## meanSD plot ##
#################

# get list of cell lines in all platforms
intersect_cells <- Reduce(intersect, lapply(names(mae), function(x){colnames(assay(mae[[x]]))})) 

mae_subset <- mae[, intersect_cells]
mae_subset@ExperimentList

# rna
p1 <- meanSdPlot(assay(mae_subset@ExperimentList$RNA), xlab  = "Rank(Mean)", ylab  = "SD") 
p1 <- p1$gg + ggtitle("Transcriptomics")  + scale_y_continuous(limits = c(0, 4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

# atac
p2 <- meanSdPlot(assay(mae_subset@ExperimentList$ATAC), xlab  = "Rank(Mean)", ylab  = "SD") 
p2 <- p2$gg + ggtitle("Epigenomics")  + scale_y_continuous(limits = c(0, 4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

# EP
p3 <- meanSdPlot(assay(mae_subset@ExperimentList$EP), xlab  = "Rank(Mean)", ylab  = "SD") 
p3 <- p3$gg + ggtitle("Expression Proteomics")  + scale_y_continuous(limits = c(0, 4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

# EP
p4 <- meanSdPlot(assay(mae_subset@ExperimentList$PP), xlab  = "Rank(Mean)", ylab  = "SD") 
p4 <- p4$gg + ggtitle("Phosphoproteomics")  + scale_y_continuous(limits = c(0, 4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

# EP
p5 <- meanSdPlot(assay(mae_subset@ExperimentList$Met), xlab  = "Rank(Mean)", ylab  = "SD") 
p5 <- p5$gg + ggtitle("Metabolomics")  + scale_y_continuous(limits = c(0, 4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/meanSD.jpeg", width = 1000, height=700)
ggarrange(p1, p2, p3, p4, p5, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()



