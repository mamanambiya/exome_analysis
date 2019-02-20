#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(plyr)
library(ggplot2)
library(reshape2)

# setwd("~/Project_data/exome_aibst/data/PGX_ONLY/FST/AIBST")
superpop.plot.colours <- c("#009E73","#CC79A7", "#56B4E9","#E69F00","#000000", "#CC0000", "#006600", "#669999", "#00CCCC",
"#660099", "#CC0066", "#FF9999", "#FF9900")

###
###    Step 14.2: Clinical variants analysis.
###
# ANALYSE CLINICAL VARIANTS
ClinPerSample_012_file <- "${ClinPerSample_012_file}"
ClinPerSample_012_indv_file <- "${ClinPerSample_012_indv_file}"
clinSampleTable_csv_file <- "${clinSampleTable_csv_file}"
clinPopTable_csv_file <- "${clinPopTable_csv_file}"
pgxClinPerSamplePerPopulation_tiff_file <- "${pgxClinPerSamplePerPopulation_tiff_file}"
dataset_sample_file <- "${dataset_sample_file}"
dataset_allLofPop_file <- "${dataset_allLofPop_file}"


# ADD POPULATION DATA  from LOF analysis
allLofPop <- read.table(dataset_allLofPop_file, sep = ',', header = T)



snpScoresClin <- read.table(ClinPerSample_012_file)
# Remove the first column which is the index of sample
snpScoresClin <- snpScoresClin[,-1]

# Flip alleles for those where the non-reference allele is the major allele
# Recode as "a" and "b" first to avoid ambiguity
for (column in 1:ncol(snpScoresClin)) {
    if (colSums(snpScoresClin[column]) > nrow(snpScoresClin)) {
        snpScoresClin[,column][snpScoresClin[,column]==2] <- "a"
        snpScoresClin[,column][snpScoresClin[,column]==0] <- "b"
        snpScoresClin[,column][snpScoresClin[,column]=="a"] <- 0
        snpScoresClin[,column][snpScoresClin[,column]=="b"] <- 2
        snpScoresClin[,column] <- as.numeric(snpScoresClin[,column])
    }
}

snpScoresClin\$totalClin <- apply(snpScoresClin,1,sum)
# Upload individuals
individualsClin <- read.table(ClinPerSample_012_indv_file)
# Make a data frame with individuals and total clinical PGx vaiants they carry
allClin <- cbind(individualsClin, snpScoresClin\$totalClin)
colnames(allClin) <- c("sample", "totalClin")
allLofClinPop <- merge(allLofPop, allClin, by="sample")
## Number of clinical variants per sample
allLofClinPop\$totalClin <- abs(allLofClinPop\$totalClin)
write.csv(allLofClinPop, clinSampleTable_csv_file)
# summary(allLofClinPop\$totalClin)
clinPopTable <- table(allLofClinPop\$totalClin,allLofClinPop\$pop)
# clinPopTable
### Fnumber of sample having y-axis clinical variants per populations
write.csv(clinPopTable, clinPopTable_csv_file)

# PLOT CLINICAL VARIANTS / SAMPLE
# Violin and jitter plot: Number of clinically-relevant variants/sample/population
clinPerSamplePlot <- ggplot(allLofClinPop, aes(x=pop, y=totalClin)) +
    geom_violin(aes(colour=dataset)) + theme_bw() +
    theme(legend.title=element_blank(),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold")) +
    xlab("Population") + ylab("Count") +
    scale_colour_manual(values=superpop.plot.colours) +
    geom_jitter(alpha=0.5, aes(color=dataset),
    position = position_jitter(width = 0.1))
# clinPerSamplePlot
ggsave(pgxClinPerSamplePerPopulation_tiff_file,  height=5.36, width=7.58, units='in', dpi=120)

## TODO
# plot_grid(clinFreqPlot, clinPerSamplePlot, labels = c("A", "B"), align = "v", nrow=2)
# ggsave("pgxClinicalGrid.tiff", height=5.36*2, width=7.58, units='in', dpi=120)

# Compare carriers to non-carriers of LoF variants
allLofClinPop\$LofCarry <- with(allLofClinPop, ifelse(totalLof==0, "No","Yes"))
allLofClinPop\$ClinCarry <- with(allLofClinPop, ifelse(totalClin==0, "No","Yes"))

tableLoFcarry <- table(allLofClinPop\$LofCarry, allLofClinPop\$pop)
# prop.table(tableLoFcarry,2)
# table(allLofClinPop\$LofCarry)
