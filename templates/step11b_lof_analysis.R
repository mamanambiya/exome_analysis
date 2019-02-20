#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(ggplot2)
library(reshape2)
library(plyr)

###
### Step 11b: Analysis of LoF
###

# ANALYSE LoF VARIANTS
# Read in the data from vcftools --012 for all individuals
# "Genotypes are represented as 0, 1 and 2,
# where the numbers represent that number of non-reference alleles"
# Missing genotypes are represented by -1

#setwd("/Users/mamana/Project_Data/exome_aibst/data/PGX_ONLY/LOF/AIBST")
snpScoresLof <- read.table("${dataset_LofPerSample_012}")
# Remove the first column which is the index of sample
snpScoresLof <- snpScoresLof[,-1]
snpScoresLof\$totalLof <- apply(snpScoresLof,1,sum)

# Read in the individuals sheet which represents
# This file represents the individuals who were analysed
individualsLof <- read.table("${dataset_LofPerSample_012_indv}")

# Make a data frame with individuals and total LoF they carry
allLof <- cbind(individualsLof, snpScoresLof\$totalLof)
colnames(allLof) <- c("sample", "totalLof")

# ADD POPULATION DATA to the DATAFRAME
samples <- read.table("${dataset_sample}",header = F)
colnames(samples) <- c("sample", "pop", "supergroup","sex", "dataset")
# Make data frame with unique populations and their superpops
popDetail <- unique(samples[,2:3])
colnames(popDetail) <- c("pop", "supergroup")
popDetail <- popDetail[with(popDetail, order(popDetail\$pop)), ]
populations <- popDetail\$pop

# allLofPop <- merge(samples\$V1, allLof, by="sample")
allLofPop <- merge(samples, allLof, by.x="sample", by.y = 'sample')
write.csv(allLofPop, "${dataset_allLofPop}")

# MAKE PLOTS in GGPLOT
# Order for plotting in ggplot
#allLofPop\$pop <- factor(allLofPop\$pop,levels=popDetail\$pop[order(popDetail\$super_pop)])
allLofPop\$pop <- factor(allLofPop\$pop, levels=popDetail\$pop[order(popDetail\$supergroup)])


# Number of individuals carrying at least one LoF variant per pop
lofSuper <- table(allLofPop\$totalLof, allLofPop\$pop)
prop_ <- prop.table(lofSuper,2)

# TODO nto needed as dealing with only one group
# Test for differences in means using Kruskal-Wallis rank sum test
#kruskal.test(totalLof ~ super_pop, data = allLofPop)
kruskal_<- kruskal.test(totalLof ~ pop, data = allLofPop)


# Plot number of LoF variants per gene and their frequency
lofConsequences <- read.table("${dataset_LoFConsequencesIDacGENE}")
colnames(lofConsequences) <- c("UniqueID", "POS", "REF", "ALT", "AC", "KG_AF", "KG_AFR_AF", "KG_EUR_AF", "KG_AMR_AF", "KG_EAS_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_FIN_AF", "ExAC_AF", "ExAC_AFR_AF", "AGVP_AF", "SAHGP_AF", "TRYPANOGEN_AF", "Gene")


##  compute sum per gene
## We take the mean frequency for a gene
lofConsequences[6:18] <- lapply(lofConsequences[6:18], as.character)
lofConsequences[6:18] <- lapply(lofConsequences[6:18], as.numeric)
lofConsequences[6:18][is.na(lofConsequences[6:18])] <- 0

## Using mean
lofACperGene_mean <- ddply(lofConsequences, c("Gene"), summarise, N=sum(as.numeric(as.character(AC))), LOF=length(AC), KG_AF=mean(as.numeric(as.character(KG_AF))), KG_AFR_AF=mean(as.numeric(as.character(KG_AFR_AF))), KG_EUR_AF=mean(as.numeric(as.character(KG_EUR_AF))),  KG_AMR_AF=mean(as.numeric(as.character(KG_AMR_AF))), KG_EAS_AF=mean(as.numeric(as.character(KG_EAS_AF))), gnomAD_AF=mean(as.numeric(as.character(gnomAD_AF))), gnomAD_AFR_AF=mean(as.numeric(as.character(gnomAD_AFR_AF))), gnomAD_EUR_AF=mean(as.numeric(as.character(gnomAD_FIN_AF))), ExAC_AF=mean(as.numeric(as.character(ExAC_AF))), ExAC_AFR_AF=mean(as.numeric(as.character(ExAC_AFR_AF))), AGVP_AF=mean(as.numeric(as.character(AGVP_AF))), SAHGP_AF=mean(as.numeric(as.character(SAHGP_AF))), TRYPANOGEN_AF=mean(as.numeric(as.character(TRYPANOGEN_AF))))
lofACperGene_mean\$GeneVariants <- with(lofACperGene_mean, ifelse(lofACperGene_mean\$LOF>1 ,paste(lofACperGene_mean\$Gene, " (", lofACperGene_mean\$LOF, " variants)", sep=""),
paste(lofACperGene_mean\$Gene, " (", lofACperGene_mean\$LOF, " variant)", sep="")))

# Reorder factors for plotting with ggplot
lofACperGene_mean\$GeneVariants <- factor(lofACperGene_mean\$GeneVariants, levels=lofACperGene_mean[order(lofACperGene_mean\$N),"GeneVariants"])

# Bubble plot of AC/no. of LoF variants in ggplot
lofACplot <- ggplot(lofACperGene_mean, aes(x = GeneVariants, y = N) ) +
    geom_point(shape=19, aes(size=LOF)) +
    theme_bw() + xlab("Gene") + ylab("Combined allele count") +
    theme(legend.position="none", axis.text.x=element_text(angle=90,
    hjust = 1, vjust = 0.5, size=12, face=c("bold")))
#lofACplot
ggsave("${dataset_pgxLoFPerGeneCount}", height=10, width=15.58, units='in', dpi=120)

# Allele frequencies of LoF variants in different super populations
# Melt and reorder for plotting
mLofACperGene_mean <- melt(lofACperGene_mean[,c(2,4:17)], id.vars=c("GeneVariants","N"))
mLofACperGene_mean\$variable <- gsub("_AF","", mLofACperGene_mean\$variable)
#Plot combined allele frequency
lofAFplot <- ggplot(mLofACperGene_mean, aes(x = GeneVariants, y = value)) +
    geom_point(shape=19, aes(colour=variable)) +
    theme_bw() + xlab("Gene") + ylab("Combined allele frequency") +
    theme(legend.title=element_blank(), legend.position="bottom",
    axis.text.x=element_text(angle=90,
    hjust = 1, vjust = 0.5, size=12, face=c("bold")))
#lofAFplot
ggsave("${dataset_pgxLoFPerCombinedAF}", height=10, width=15.58, units='in', dpi=120)
# SLC19A1 in AMR is driven mainly by PEL individuals (rs200964006)

#plot_grid(lofACplot, lofAFplot, labels = c("A", "B"), align = "v", nrow=2)
