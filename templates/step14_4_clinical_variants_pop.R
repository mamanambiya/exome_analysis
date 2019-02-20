#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(plyr)
library(ggplot2)
library(reshape2)

# setwd("~/Project_data/exome_aibst/data/PGX_ONLY/FST/AIBST")
superpop.plot.colours <- c("#3A8B8B","#7CFC00","#CD9B1D","#7FFFD4","#FF8C00","#6495ED","#008B8B","#FF0000",
"#008B45","#8B4500","#8B7355","#CDCD00","#98FB98","#8B4513","#FFD700","#00FFFF",
"#98F5FF","#228B22","#548B54","#8B3A62","#CD6090","#CDC673","#FF69B4","#FFEC8B",
"#CD6600","#2E8B57")


###
   ## Step 14.4.
###
## ANALYSE THE CLINICALLY RELEVANT VARIANTS
# rsIDGene <- read.table("${pgxClinicalLevel1}", header = T, sep = "\\t")
# colnames(rsIDGene) <- c("UniqueID", "rsID", "GENE")
rsIDGene <- read.table("${pgxClinicalLevel1}", header = T, sep = "\\t")
rsIDGene <- rsIDGene[,1:7]
colnames(rsIDGene) <- c("CHROM", "POS", "rsID", "REF", "ALT", "AC", "GENE")

datasetAnnotated <- read.table("${datasetAnnotated}", header = T, sep = ",")
clinDataset <- merge(datasetAnnotated, rsIDGene, by="rsID")

samples <- read.table("${dataset_sample_file}", header = F)
colnames(samples) <- c("sample", "pop", "supergroup","sex", "dataset")
# Make data frame with unique populations and their superpops
popDetail <- unique(samples[,2:3])
# colnames(popDetail) <- c("pop", "supergroup")
popDetail <- popDetail[with(popDetail, order(popDetail\$pop)), ]
populations <- popDetail\$pop
myPop <- append(as.vector(popDetail\$pop), "KG_AFR")
myPop <- append(myPop, "gnomAD_AFR")
# myPop <- append(myPop, "KG")
myPop <- append(myPop, "AGVP")

# Select all columns with ref allele frequency
RefAlleleFreqClin <- clinDataset[,c(grep('*_Alt_FREQ',names(clinDataset)))]
RefAlleleFreqClin\$rsID <- clinDataset\$rsID
RefAlleleFreqClin\$Gene <- clinDataset\$GENE
RefAlleleFreqClin <- RefAlleleFreqClin[order(RefAlleleFreqClin\$Gene),]
write.csv(RefAlleleFreqClin, "${tsv_file}")


# Plot allele frequencies
# Make tall data frame for plotting in ggplot
mRefAlleleFreqClin <- melt(RefAlleleFreqClin, id=c("Gene", "rsID"))
mRefAlleleFreqClin[mRefAlleleFreqClin == '.']<- 0
# mRefAlleleFreqClin\$Pop <- rep(myPop, each=nrow(RefAlleleFreqClin))
mRefAlleleFreqClin\$Pop <- gsub("_Alt_FREQ", "",mRefAlleleFreqClin\$variable)
# mRefAlleleFreqClin\$Super <- rep(popDetail\$supergroup[order(myPop)],each=nrow(RefAlleleFreqClin))
mRefAlleleFreqClin\$rsID <- sapply(strsplit(as.character(RefAlleleFreqClin\$rsID), ";"), `[`, 1)
mRefAlleleFreqClin\$GeneRs <- paste(mRefAlleleFreqClin\$Gene, mRefAlleleFreqClin\$rsID, sep=" ")
mRefAlleleFreqClin\$superPop <- mRefAlleleFreqClin\$Pop
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='KNK', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='KNL', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='KNM', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='KNP', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='NGH', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='NGI', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='NGY', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='SAV', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='TZA', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='TZB', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='ZWD', 'AIBST_POP', mRefAlleleFreqClin\$superPop)
mRefAlleleFreqClin\$superPop <- ifelse(mRefAlleleFreqClin\$Pop=='ZWS', 'AIBST_POP', mRefAlleleFreqClin\$superPop)

# Plot allele frequency stripchart
# Allele frequencies of pharmacogenomic variants with high levels of clinical evidence
clinFreqPlot <- ggplot(mRefAlleleFreqClin, aes(x = GeneRs, y = as.numeric(value))) +
geom_point(aes(colour = Pop, shape=superPop)) +
  theme_bw() + xlab("Variant") +
  scale_colour_manual(values=superpop.plot.colours) +
  ylab("Frequency") +
scale_y_continuous(breaks=seq(0, 0.5, 0.1)) +
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1,
vjust = 0.5, size=12))
clinFreqPlot
ggsave("${tiff_file}",
height=6.20, width=13, units='in', dpi=120)
