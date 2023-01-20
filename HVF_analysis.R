########################################
##   Investigation of HVF selection   ##
########################################

###############
## libraries ##
###############

library(tidyverse)
source("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/MODI_helper_functions.R")

###############
## load data ##
###############

mae <- readRDS("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/mae_all_2022.08.05_naomit.rds")
intersect_cells <- Reduce(intersect, lapply(names(mae), function(x){colnames(assay(mae[[x]]))})) 
mae_subset <- mae[, intersect_cells]

mae2 <- readRDS("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/mae_all_2022.09.01.rds")

###########################
## feature variance plot ##
###########################

d1 <- data.frame("Variance" = sort(matrixStats::rowVars(assay(mae_subset[["RNA"]])), decreasing=TRUE)) %>%
  rownames_to_column("Rank") %>%
  mutate(Rank=as.numeric(Rank))
d1_x <- d1 %>%
  filter(Rank==5000) %>%
  pull(Variance)

p1 <- ggplot(d1) +
  geom_point(aes(x=Rank, y=Variance, colour = Rank<=5000)) +
  scale_colour_manual(name = 'Feature Included', values = setNames(c('#2266ac','#b2182b'),c(T, F))) +
  geom_vline(xintercept = 5000, linetype="dashed") +
  geom_hline(yintercept=d1_x, linetype="dashed") +
  ggtitle("Transcriptomics") + 
  scale_y_continuous(limits = c(0, 16), breaks=seq(0,16,2)) + 
  scale_x_continuous(breaks=seq(0,30000,5000)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

d2 <- data.frame("Variance" = sort(matrixStats::rowVars(assay(mae_subset[["ATAC"]])), decreasing=TRUE)) %>%
  rownames_to_column("Rank") %>%
  mutate(Rank=as.numeric(Rank))
d2_x <- d2 %>%
  filter(Rank==5000) %>%
  pull(Variance)

p2 <- ggplot(d2) +
  geom_point(aes(x=Rank, y=Variance, colour = Rank<=5000)) +
  scale_colour_manual(name = 'Feature Included', values = setNames(c('#2266ac','#b2182b'),c(T, F))) +
  geom_vline(xintercept = 5000, linetype="dashed") +
  geom_hline(yintercept=d2_x, linetype="dashed") +
  ggtitle("Epigenomics") + 
  scale_y_continuous(limits = c(0, 16), breaks=seq(0,16,2)) + 
  scale_x_continuous(breaks=seq(0,33000,5000)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

d3 <- data.frame("Variance" = sort(matrixStats::rowVars(assay(mae_subset[["EP"]])), decreasing=TRUE)) %>%
  rownames_to_column("Rank") %>%
  mutate(Rank=as.numeric(Rank))
d3_x <- d3 %>%
  filter(Rank==2000) %>%
  pull(Variance)

p3 <- ggplot(d3) +
  geom_point(aes(x=Rank, y=Variance, colour = Rank<=2000)) +
  scale_colour_manual(name = 'Feature Included', values = setNames(c('#2266ac','#b2182b'),c(T, F))) +
  geom_vline(xintercept = 2000, linetype="dashed") +
  geom_hline(yintercept=d3_x, linetype="dashed") +
  ggtitle("Expression Proteomics") + 
  scale_y_continuous(limits = c(0, 16), breaks=seq(0,16,2)) + 
  scale_x_continuous(breaks=seq(0,6500,2000)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

d4 <- data.frame("Variance" = sort(matrixStats::rowVars(assay(mae_subset[["PP"]])), decreasing=TRUE)) %>%
  rownames_to_column("Rank") %>%
  mutate(Rank=as.numeric(Rank))
d4_x <- d4 %>%
  filter(Rank==2000) %>%
  pull(Variance)

p4 <- ggplot(d4) +
  geom_point(aes(x=Rank, y=Variance, colour = Rank<=2000)) +
  scale_colour_manual(name = 'Feature Included', values = setNames(c('#2266ac','#b2182b'),c(T, F))) +
  geom_vline(xintercept = 2000, linetype="dashed") +
  geom_hline(yintercept=d4_x, linetype="dashed") +
  ggtitle("Phosphoproteomics") + 
  scale_y_continuous(limits = c(0, 16), breaks=seq(0,16,2)) + 
  scale_x_continuous(breaks=seq(0,5000,2000)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

d5 <- data.frame("Variance" = sort(matrixStats::rowVars(assay(mae_subset[["Met"]])), decreasing=TRUE)) %>%
  rownames_to_column("Rank") %>%
  mutate(Rank=as.numeric(Rank))
d5_x <- d5 %>%
  filter(Rank==1000) %>%
  pull(Variance)

p5 <- ggplot(d5) +
  geom_point(aes(x=Rank, y=Variance, colour = Rank<=1000)) +
  scale_colour_manual(name = 'Feature Included', values = setNames(c('#2266ac','#b2182b'),c(T, F))) +
  geom_vline(xintercept = 1000, linetype="dashed") +
  geom_hline(yintercept=d5_x, linetype="dashed") +
  ggtitle("Metabolomics") + 
  scale_y_continuous(limits = c(0, 16), breaks=seq(0,16,2)) + 
  scale_x_continuous(breaks=seq(0,1500,500)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/feature_variability.jpeg", width = 1000, height=700)
ggarrange(p1, p2, p3, p4, p5, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

##################
## HVF overlaps ##
##################

## get lists of HVF for each omics

# filtered
top_varFeats_filtered <- list("ATAC" = top_varFeats_view(mae_subset, "ATAC", NA_omit=TRUE, 5000),
                     "EP" = top_varFeats_view(mae_subset, "EP", NA_omit=TRUE, 2000),
                     "Met" = top_varFeats_view(mae_subset, "Met", NA_omit=TRUE, 1000),
                     "PP" = top_varFeats_view(mae_subset, "PP", NA_omit=TRUE, 2000),
                     "RNA" = top_varFeats_view(mae_subset, "RNA", NA_omit=TRUE, 5000))

# unfiltered
top_varFeats_unfiltered <- list("ATAC" = top_varFeats_view(mae2, "ATAC", NA_omit=FALSE, 5000),
                     "EP" = top_varFeats_view(mae2, "EP", NA_omit=FALSE, 2000),
                     "Met" = top_varFeats_view(mae2, "Met", NA_omit=FALSE, 1000),
                     "PP" = top_varFeats_view(mae2, "PP", NA_omit=FALSE, 2000),
                     "RNA" = top_varFeats_view(mae2, "RNA", NA_omit=FALSE, 5000))

# plot venn diagrams
# library
library(VennDiagram)

# Chart
venn.diagram(
  x = list(top_varFeats_unfiltered[["RNA"]] %>% unlist() , 
           top_varFeats_filtered[["RNA"]] %>% unlist()),
  category.names = c("Unfiltered" , "Filtered"),
  filename = 'RNA.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.035, 0.035)
)

venn.diagram(
  x = list(top_varFeats_unfiltered[["ATAC"]] %>% unlist() , 
           top_varFeats_filtered[["ATAC"]] %>% unlist()),
  category.names = c("Unfiltered" , "Filtered"),
  filename = 'ATAC.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.035, 0.035)
)
venn.diagram(
  x = list(top_varFeats_unfiltered[["EP"]] %>% unlist() , 
           top_varFeats_filtered[["EP"]] %>% unlist()),
  category.names = c("Unfiltered" , "Filtered"),
  filename = 'EP.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.035, 0.035)
)
venn.diagram(
  x = list(top_varFeats_unfiltered[["PP"]] %>% unlist() , 
           top_varFeats_filtered[["PP"]] %>% unlist()),
  category.names = c("Unfiltered" , "Filtered"),
  filename = 'PP.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.035, 0.035)
)
venn.diagram(
  x = list(top_varFeats_unfiltered[["Met"]] %>% unlist() , 
           top_varFeats_filtered[["Met"]] %>% unlist()),
  category.names = c("Unfiltered" , "Filtered"),
  filename = 'Met.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.035, 0.035)
)
# look at difference in min var 

RNA <- assay(mae_subset[["RNA"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
RNA <- min(head(sort(RNA, decreasing = TRUE), 5000))

ATAC <- assay(mae_subset[["ATAC"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
ATAC <- min(head(sort(ATAC, decreasing = TRUE), 5000))

EP <- assay(mae_subset[["EP"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
EP <- min(head(sort(EP, decreasing = TRUE), 2000))

PP <- assay(mae_subset[["PP"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
PP <- min(head(sort(PP, decreasing = TRUE), 2000))

Met <- assay(mae_subset[["Met"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
Met <- min(head(sort(Met, decreasing = TRUE), 1000))

RNA2 <- assay(mae2[["RNA"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
RNA2 <- min(head(sort(RNA2, decreasing = TRUE), 5000))

ATAC2 <- assay(mae2[["ATAC"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
ATAC2 <- min(head(sort(ATAC2, decreasing = TRUE), 5000))

EP2 <- assay(mae2[["EP"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
EP2 <- min(head(sort(EP2, decreasing = TRUE), 2000))

PP2 <- assay(mae2[["PP"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
PP2 <- min(head(sort(PP2, decreasing = TRUE), 2000))

Met2 <- assay(mae2[["Met"]]) %>% 
  as.matrix() %>% 
  matrixStats::rowVars(.)
Met2 <- min(head(sort(Met2, decreasing = TRUE), 1000))


