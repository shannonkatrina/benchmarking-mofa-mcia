####################################################################
## Filter only promoter peaks from ATAC-seq data using ChIPseeker ##
####################################################################

###############
## libraries ##
###############

library(ChIPseeker)
library(tidyverse)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
data_dir <- "/home/skb45807/MODI/data/47_cancer_cells/"

###############
## load data ##
###############

atac_data <- read.table(paste0(data_dir,"ATAC-seq/Combined_Filtered_BrCa_LuCa_OvCa_RawReads8_Deseq2LibSize_VST_Normalized_Median_Peak_Counts_LongFormat.tsv"), header = T)
all_peaks <- atac_data %>%
  dplyr::select(PeakID) %>%
  separate(PeakID,sep="_",into=c("chr","start","end")) %>%
  mutate(start=as.numeric(start),
         end=as.numeric(end)) %>%
  na.omit()

# convert tsv to bed files
write.table(all_peaks, paste0(data_dir,"ATAC-seq/peaks.bed"),  sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

####################
## annotate peaks ##
####################

all_peaks_GRanges <- readPeakFile(paste0(data_dir,"ATAC-seq/peaks.bed")) # 3 column bed file with Chr, St, End
peakAnno <- annotatePeak(all_peaks_GRanges, tssRegion=c(-200, 200), addFlankGeneInfo=TRUE, flankDistance=200,
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno) # ~8% peaks in promoter region

############################
## extract promoter peaks ##
############################

peakAnno@anno$annotation #peak annotation
peakAnno@anno$ENSEMBL #peak id

data.frame(ensemblID=peakAnno@anno$ENSEMBL,
           annotation=peakAnno@anno$SYMBOL)

peaksGranges_Anno <- as.GRanges(peakAnno)
promoter_partIDs <- as.data.frame(peaksGranges_Anno[peaksGranges_Anno$annotation == "Promoter"]) %>%
  mutate(start=start-1, #some reason all the start ranges are +1 vs starts in all peaks
    ID = paste(seqnames,start,end,sep="_")) %>%
  pull(ID)

promoter_IDs <- atac_data %>%
  mutate(partID = sub("_peak.*", "", PeakID)) %>%
  dplyr::filter(partID %in% promoter_partIDs) %>%
  pull(PeakID)

# save list of promoter peak IDs 
write.table(promoter_IDs, paste0(data_dir,"ATAC-seq/promoter_PeakIDs.txt"), sep="\t", col.names = FALSE, row.names = FALSE)

########################################

#########################
## downstream analysis ##
#########################

## filter 4 peaks of interest
peaks_of_interest <- as.data.frame(peaksGranges_Anno[peaksGranges_Anno$annotation == "Promoter"]) %>%
  mutate(start=start-1, #some reason all the start ranges are +1 vs starts in all peaks
         ID = paste(seqnames,start,end,sep="_")) %>%
  filter(ID %in% c('chr8_101368515_101369549','chr14_23153980_23154829','chr3_191360615_191360975','chr18_35345217_35345972')) %>%
  dplyr::select(ID, SYMBOL)


