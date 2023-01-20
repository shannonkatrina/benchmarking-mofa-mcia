#############################################
##   Comparison of MOFA and MCIA outputs   ##
##             R version 4.0.2             ##
#############################################

###############
## libraries ##
###############

library(tidyverse)
library(MOFA2)
library(rlist)
library(reticulate)
library(omicade4)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
source("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/MODI_helper_functions.R")
use_condaenv("/Users/shannonbrown/opt/anaconda3/envs/mofa2/", required=TRUE) 
py_config()

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

###############
## load data ##
###############

MOFAobject.trained <- readRDS("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/MOFAtrained.rds")

mcoin <- readRDS("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/MCOIN.rds")

metadata <- read.csv("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/new_cell_metadata_2022.08.03.csv", na.strings=c("","NA")) 

mae <- readRDS("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/mae_all_2022.08.05_naomit.rds")
intersect_cells <- Reduce(intersect, lapply(names(mae), function(x){colnames(assay(mae[[x]]))})) 

drug_response<- read.csv("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/PANCANCER_IC_Fri Aug  5 16_49_38 2022.csv") 
drug_response_subset <- drug_response %>%
  mutate(stripped_cell_line_name = toupper(str_replace_all(Cell.Line.Name, "[[:punct:]]", "")),
         drug = paste0(Drug.Name,"_",Drug.ID,"_",Max.Conc))  %>%
  filter(stripped_cell_line_name %in% intersect_cells) %>%
  na.omit()

############################
## extracts factor values ##
############################

mofa_eigs <- as.data.frame(get_factors(MOFAobject.trained)[[1]]) %>%
  rownames_to_column("stripped_cell_line_name")

mcoin_eigs <- as.data.frame(mcoin$mcoa$SynVar) %>%  #synthetic score (Factor) for each cell line 
  rownames_to_column("stripped_cell_line_name") 
colnames(mcoin_eigs) <- gsub("SynVar","Factor",colnames(mcoin_eigs))


#######################################
## variance explained visualisations ##
#######################################

## variance decomposition ##

# MOFA
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/var_decomp_2.jpeg", width = 600, height=400)
plot_variance_explained(MOFAobject.trained, x="view", y="factor")
dev.off()

# MCIA
factors <- mcoin$mcoa$cov2 %>%
  rownames_to_column(var="Omics") %>%
  pivot_longer(cols=-Omics, names_to = "Factor", values_to = "Var. (%)") %>%
  mutate("Var. (%)"=round(`Var. (%)`*100, digits=2),
         Factor=gsub("cov2","Factor",Factor))

factors$Factor <- factor(factors$Factor, levels = paste0(rep("Factor",15),1:15))
factors$Omics <- factor(factors$Omics, levels = c("ATAC","EP","PP","RNA"))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/var_decomp_2.jpeg", width = 600, height=400)
ggplot(factors, aes(x=Omics, y=Factor, fill=`Var. (%)`)) +
  geom_tile(colour="black", stat="identity") +
  scale_fill_gradient(low = "white", high = "#1f1095", limits=c(0,max(factors$`Var. (%)`))) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, color="black", family="Arial"),
        axis.text.y = element_text(size=12, color="black", family="Arial")) +
  xlab("") +
  ylab("")
dev.off()

## total variance explained ##

# MOFA
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/var_decomp_omics_2.jpeg", width = 600, height=400)
plot_variance_explained(MOFAobject.trained, x="group", y="factor", plot_total = T)[[2]]
dev.off()

# MCIA
omics_total <- factors %>%
  group_by(Omics) %>%
  summarise(`Variance explained (%)` = sum(`Var. (%)`))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/var_decomp_omics_2.jpeg", width = 600, height=400)
ggplot(omics_total, aes(x=Omics, y=`Variance explained (%)`)) +
  geom_bar(stat="identity", fill="#02698b", colour="black", width = 0.75) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=12, color="black", family="Arial"),
        axis.text.y = element_text(size=12, color="black", family="Arial"),
        axis.title.y = element_text(size=12, color="black", family="Arial"),
        axis.line = element_line(size=rel(1.0), color="black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=rel(1.0))) +
  xlab("") +
  scale_x_discrete(limits=omics_total$Omics)
dev.off()

## intra-method factor correlation ##

# MOFA
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/corr_2.jpeg", width = 500, height=500)
plot_factor_cor(MOFAobject.trained)
dev.off()

# MCIA
M <- cor(mcoin_eigs[,-1])
jpeg("/home/skb45807/MODI/visualisations/47_cancer_cells/MCIA/corr_2.jpeg", width = 500, height=500)
corrplot(M)
dev.off()

#########################
## factor-factor plots ## 
##     cancer type     ##
#########################

## factor 1 vs factor 2 ##

# MOFA
# add metadata to model
model_metadata <- metadata %>%
  rename(sample = stripped_cell_line_name) %>%
  filter(sample %in% MOFAobject.trained@samples_metadata$sample)
samples_metadata(MOFAobject.trained) <- model_metadata ### can't add unless all metadata there
head(MOFAobject.trained@samples_metadata,3) # additional fields added

# plot
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/LF1vsLF2_cancer.jpeg", width = 600, height=500)
plot_factors(MOFAobject.trained, 
  factors = 1:2,
  color_by = "primary_disease"
) +
  geom_point(size=3)
dev.off()

# MCIA
mcoin_eigs_anno <- mcoin_eigs %>%
  left_join(metadata %>% select(stripped_cell_line_name, primary_disease, Subtype), by="stripped_cell_line_name")
  
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/LF1vsLF2_cancer.jpeg", width = 600, height=500)
ggplot(mcoin_eigs_anno, aes(x=Factor1, y=Factor2)) +
  geom_point(aes(colour=primary_disease), size=3) +
  geom_point(shape=1, size=3, colour="black") +
  theme(axis.text.x = element_text(size=12, color="black", family="Arial"),
        axis.text.y = element_text(size=12, color="black", family="Arial"),
        axis.title.y = element_text(size=12, color="black", family="Arial"),
        axis.title.x = element_text(size=12, color="black", family="Arial"),
        axis.line = element_line(size=rel(1.0), color="black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=rel(1.0)))
dev.off()

## all combinations for supplement ##

# MOFA
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/LFpairwise_cancer.jpeg", width = 900, height=800)
plot_factors(MOFAobject.trained, 
  factors = 1:10,
  color_by = "primary_disease"
)
dev.off()

# MCIA
library(GGally)
jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/LFpairwise_cancer.jpeg", width = 900, height=800)
ggpairs(mcoin_eigs_anno,
        columns=2:11,aes(color=primary_disease),
        upper=list(continuous = "points"),
        legend = 2) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


#####################################
## compare corr of method factors ##
#####################################

all_eigs_wide <- left_join(mofa_eigs,mcoin_eigs, by="stripped_cell_line_name")
M = cor(all_eigs_wide[,-1])

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/mcia_mofa_corrplot.jpeg", width = 1100, height=1000)
corrplot(M) # by default, method = 'circle'
dev.off()

# plot specific factor comparisons
mcoin_eigs_subset <- mcoin_eigs %>%
  select(stripped_cell_line_name, mcia_1=Factor1, mcia_5=Factor5)
mofa_eigs_subset <- mofa_eigs %>%
  select(stripped_cell_line_name, mofa_1=Factor1, mofa_6=Factor6)

all_eigs <- left_join(mcoin_eigs_subset,mofa_eigs_subset, by="stripped_cell_line_name")

p1 <- ggplot(all_eigs, aes(x=mofa_1, y=mcia_1)) +
  geom_point() +
  geom_smooth(method=lm, se=TRUE) +
  stat_cor(method = "pearson") +
  ylab("MCIA Factor 1") +
  xlab("MOFA Factor 1") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

p2 <- ggplot(all_eigs, aes(x=mofa_6, y=mcia_5)) +
  geom_point() +
  geom_smooth(method=lm, se=TRUE) +
  stat_cor(method = "pearson") +
  ylab("MCIA Factor 5") +
  xlab("MOFA Factor 6") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/mofa_mcia_scatter.jpeg", width = 900, height=600)
ggarrange(p1, p2, ncol=2)
dev.off()

######################################
## drug response variance explained ##
######################################

## filter drugs that have at least 20 observations
drug_filter <- drug_response_subset %>%
  group_by(drug) %>%
  summarise(n=n()) %>%
  filter(n>=20) %>%
  pull(drug)

# functions for calculating R squared and p-value
rsq <- function(x, y){summary(lm(y~x))$r.squared}
lmp <- function (x, y) {
  modelobject <- lm(y~x)
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# colour pallete of blues for 10 factors
coul <- brewer.pal(4, "Blues") 
coul <- colorRampPalette(coul)(10) 

## MOFA ##

mofa_eigs_long_anno <- mofa_eigs %>%
  pivot_longer(cols=-stripped_cell_line_name, names_to="factor", values_to="eigenvalue") %>%
  left_join(drug_response_subset %>% select(stripped_cell_line_name, drug, AUC), by="stripped_cell_line_name")

mofa_drug_rsq <- mofa_eigs_long_anno %>% 
  select(-stripped_cell_line_name) %>%
  filter(drug %in% drug_filter) %>%
  group_by(factor, drug) %>%
  summarize(rsquared=rsq(eigenvalue, AUC),
            pval=lmp(eigenvalue, AUC)) %>%
  ungroup() %>%
  group_by(factor) %>%
  mutate(p.adj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  rename(Factor=factor) %>%
  mutate(drug_name=sub("_.*", "", drug))

mofa_drug_corr <- mofa_eigs_long_anno %>% 
  select(-stripped_cell_line_name) %>%
  filter(drug %in% drug_filter) %>%
  group_by(factor, drug) %>%
  summarize(COR = stats::cor.test(eigenvalue, AUC)$estimate,
            pval = stats::cor.test(eigenvalue, AUC)$p.value
  ) %>%
  mutate(p.adj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  rename(Factor=factor) %>%
  mutate(drug_name=sub("_.*", "", drug))

# boxplot of r-sqaured/corr for each factor with annotated outliers
mofa_drug_rsq$Factor <- factor(mofa_drug_rsq$Factor, levels = paste0(rep("Factor",15),1:15))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/mofa_rsq_boxplot.jpeg", width = 800, height=400)
ggplot(mofa_drug_rsq, aes(x=Factor, y=rsquared, fill=Factor, label=drug_name)) +
  geom_boxplot() +
  geom_text(aes(label=ifelse(rsquared>0.3,as.character(drug_name),'')),hjust=0,vjust=-1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"),
        legend.position = "none") +
  ylab("R²") +
  xlab("") +
  scale_fill_manual(values=rev(dput(as.character(coul))))
dev.off()

ggplot(mofa_drug_corr, aes(x=Factor, y=COR, fill=Factor, label=drug_name)) +
  geom_boxplot() +
  geom_text(aes(label=ifelse(COR>0.5,as.character(drug_name),'')),hjust=0,vjust=-1) +
  geom_hline(yintercept = c(-0.5,0.5), linetype="dashed") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"),
        legend.position = "none") +
  ylab("Correlation") +
  xlab("")+
  scale_fill_manual(values=rev(dput(as.character(coul))))


## MCIA ##

mcoin_eigs_long_anno <- mcoin_eigs %>%
  pivot_longer(cols=-stripped_cell_line_name, names_to="factor", values_to="eigenvalue") %>%
  left_join(drug_response_subset %>% select(stripped_cell_line_name, drug, AUC), by="stripped_cell_line_name")

mcoin_drug_rsq <- mcoin_eigs_long_anno %>% 
  select(-stripped_cell_line_name) %>%
  filter(drug %in% drug_filter) %>%
  group_by(factor, drug) %>%
  summarize(rsquared=rsq(eigenvalue, AUC),
            pval=lmp(eigenvalue, AUC)) %>%
  ungroup() %>%
  group_by(factor) %>%
  mutate(p.adj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  rename(Factor=factor) %>%
  mutate(drug_name=sub("_.*", "", drug))

# boxplot of r-sqaured/corr for each factor with annotated outliers
mcoin_drug_rsq$Factor <- factor(mcoin_drug_rsq$Factor, levels = paste0(rep("Factor",15),1:15))

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MCIA/mcia_rsq_boxplot.jpeg", width = 800, height=400)
ggplot(mcoin_drug_rsq, aes(x=Factor, y=rsquared, fill=Factor, label=drug_name)) +
  geom_boxplot() +
  geom_text(aes(label=ifelse(rsquared>0.3,as.character(drug_name),'')),hjust=0,vjust=-1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size=rel(1.0), color="black"),
        legend.position = "none") +
  ylab("R²") +
  xlab("") +
  scale_fill_manual(values=rev(dput(as.character(coul))))
dev.off()

#################################
## investigate AZD7762 in MOFA ##
#################################

## scatter plot of AUC vs factor value ##

AZD7762_eigsxdrugs <- mofa_eigs_long_anno %>%
  filter(drug == "AZD7762_1022_10" & factor=="Factor8") %>%
  mutate(AUC.z = (AUC-mean(AUC))/sd(AUC),
         upper_limit = 0.8*sd(AUC) + mean(AUC),
         lower_limit = -0.8*sd(AUC) + mean(AUC)
  ) %>%
  mutate(`Response Category` = case_when(
    AUC.z <= -0.8 ~ "Responder",
    AUC.z >= 0.8 ~ "Non-Responder",
    TRUE ~ "Intermediate"
  )) %>%
  left_join(metadata %>% select(stripped_cell_line_name,`Primary Disease`=primary_disease), by="stripped_cell_line_name")

table(AZD7762_eigsxdrugs$`Response Category`)

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/AZD7762_scatter.jpeg", width = 600, height=500)
ggplot(AZD7762_eigsxdrugs,aes(x=eigenvalue,y=AUC)) +
  geom_point(aes(colour=`Response Category`, shape=`Primary Disease`), size=3) +
  geom_smooth(method=lm, se=TRUE, color="dark grey") +
  stat_cor(method = "pearson") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line(size=rel(1.0), color="black")) +
  xlab("Eigenvalue") +
  ylab("AUC") 
dev.off()

## density plot of AUC for 38 cell lines ##

drug_response_AZD7762 <- drug_response_subset %>%
  filter(Drug.Name == "AZD7762" & stripped_cell_line_name %in% intersect_cells)

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/AZD7762_density.jpeg", width = 500, height=500)
density_plot = ggplot(drug_response_AZD7762, aes(x=AUC)) + 
  geom_density(alpha=0.3) + 
  geom_vline(aes(xintercept=mean(AUC)),linetype="dashed", color="black") +
  geom_vline(aes(xintercept=0.8*sd(AUC) + mean(AUC)),linetype="dotted", color="black") +
  geom_vline(aes(xintercept=-0.8*sd(AUC) + mean(AUC)),linetype="dotted", color="black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=12, color="black", family="Arial"),
        axis.text.y = element_text(size=12, color="black", family="Arial"),
        axis.line = element_line(size=rel(1.0), color="black"),
        legend.title = element_blank(),
        legend.position = "none") +
  ylab("Density") + 
  scale_fill_brewer(palette = "Paired")

dpb <- ggplot_build(density_plot)

x1 <- min(which(dpb$data[[1]]$x >= 0))
x2 <- max(which(dpb$data[[1]]$x <= 0.612913))
x3 <- min(which(dpb$data[[1]]$x >= 0.612913))
x4 <- max(which(dpb$data[[1]]$x <= 0.8058634))
x5 <- min(which(dpb$data[[1]]$x >= 0.8058634))
x6 <- max(which(dpb$data[[1]]$x <= 1))

density_plot +
  geom_area(data=data.frame(x=dpb$data[[1]]$x[x1:x2],
                            y=dpb$data[[1]]$y[x1:x2]),
            aes(x=x, y=y), fill="#619CFF", alpha = 0.5) +
  geom_area(data=data.frame(x=dpb$data[[1]]$x[x3:x4],
                            y=dpb$data[[1]]$y[x3:x4]),
            aes(x=x, y=y), fill="#F8766D", alpha = 0.5) +
  geom_area(data=data.frame(x=dpb$data[[1]]$x[x5:x6],
                            y=dpb$data[[1]]$y[x5:x6]),
            aes(x=x, y=y), fill="#00BA38", alpha = 0.5) +
  scale_color_manual(name='Regression Model',
                     breaks=c('Linear', 'Quadratic', 'Cubic'),
                     values=c('Cubic'='pink', 'Quadratic'='blue', 'Linear'='purple'))

dev.off()

## top epigenomic feature weights for mofa factor 8 ##

p1 <- plot_top_weights(MOFAobject.trained,
                      view = "ATAC",
                      factor = 8,
                      nfeatures = 10
)

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/factor8_ATACWeights.jpeg", width = 600, height=600)
p1
dev.off()

## top epigenomic features correlation for mofa factor 8 ##

model_metadata <- metadata %>%
  rename(sample = stripped_cell_line_name) %>%
  filter(sample %in% intersect_cells) %>%
  left_join(.,almost_sig %>% select(stripped_cell_line_name, `Response Category`), by=c("sample"="stripped_cell_line_name"))
samples_metadata(MOFAobject.trained) <- model_metadata
samples_metadata(MOFAobject.trained)

p2 <- plot_data_scatter(MOFAobject.trained, 
                        view = "ATAC",
                        color_by = "Response Category",
                        factor = 8,  
                        features = 4
) + labs(y="Peak Count")

jpeg("/Users/shannonbrown/OneDrive - GSK/Documents/mRes Course/mRes Project/MOFA/factor8_ATAC.jpeg", width = 1200, height=600)
ggarrange(p1, p2, ncol=2, legend="bottom")
dev.off()


