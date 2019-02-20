#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(plyr)
library(ggplot2)
library(reshape2)
library(Hmisc)

# setwd("~/Project_data/exome_aibst/data/PGX_ONLY/FST/AIBST")
superpop.plot.colours <- c("#009E73","#CC79A7", "#56B4E9","#E69F00","#000000", "#CC0000", "#006600", "#669999", "#00CCCC",
"#660099", "#CC0066", "#FF9999", "#FF9900")

###
###    Step 16.2: Consequences.
###

dataset_SO_terms_MAF_summary <- "${dataset_SO_terms_MAF_summary}"
dataset_pgxFunctionalClassesCounts_tiff <- "${dataset_pgxFunctionalClassesCounts_tiff}"
dataset_pgxFunctionalClassesByFrequency_tiff <- "${dataset_pgxFunctionalClassesByFrequency_tiff}"

# Analyse the functional annotations of variants
functionalSO <- read.table(dataset_SO_terms_MAF_summary, header=T)
# Drop 0 count terms
functionalSO <- functionalSO[functionalSO\$countAll!=0,]
# Tidy up terms for plotting
functionalSO\$terms <- gsub("_", " ", functionalSO\$terms)
functionalSO\$terms <- gsub(" variant", "", functionalSO\$terms)
functionalSO\$terms <- gsub("non coding", "NC", functionalSO\$terms)
functionalSO\$terms <- Hmisc::capitalize(functionalSO\$terms)
# Reorder based on count
functionalSO <- functionalSO[with(functionalSO, order(functionalSO\$countAll)), ]
functionalSO\$terms <- factor(functionalSO\$terms, level=functionalSO[order(functionalSO\$countAll),"terms"])
rownames(functionalSO) <- NULL

# Plot the total counts for each of the functional classes
functionalSOplot <- ggplot(functionalSO, aes(x = terms, y = countAll) ) +
    geom_point(shape=19) +
    theme_bw() + xlab("Consequence") +  ylab("Count") +
    annotate("text", size=3, fontface="bold", y=functionalSO\$countAll+2800,
    x=1:nrow(functionalSO), label=functionalSO\$countAll) +
    theme(legend.position="none", axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face=c("bold")))
functionalSOplot
ggsave(dataset_pgxFunctionalClassesCounts_tiff, height=5.36*0.75, width=7.58, units='in', dpi=120)

#Make Frequency dataframe
functionalSOFreq <- NULL
functionalSOFreq\$terms <- functionalSO\$terms
functionalSOFreq\$freqAll <- functionalSO\$countAll/sum(functionalSO\$countAll)
functionalSOFreq\$freqSingletons <- functionalSO\$countSingletons/sum(functionalSO\$countSingletons)
functionalSOFreq\$freqnonSing.0.01MAF <- functionalSO\$countnonSing.0.01MAF/sum(functionalSO\$countnonSing.0.01MAF) #To reduce size of data for Fisher exact test
functionalSOFreq\$freq0.01.0.05MAF <- functionalSO\$count0.01.0.05MAF/sum(functionalSO\$count0.01.0.05MAF)
functionalSOFreq\$freqGreater0.05 <- functionalSO\$countGreater0.05/sum(functionalSO\$countGreater0.05)
functionalSOFreq <- as.data.frame(functionalSOFreq)
# Add spaces in names for better plotting
# colnames(functionalSOFreq) <- c("terms", "freqAll", "Singletons", "Singletons - MAF 0.01", "MAF 0.01 - 0.05","MAF greater than 0.05")
colnames(functionalSOFreq) <- c("terms", "freqAll", "Singletons", "Singletons - MAF 0.01", "MAF 0.01 - 0.05", "MAF greater than 0.05")

# Test for differences between the different counts in the different frequency classes
# Make a dataframe to put significance test results in a set a count for indexing
functionalSignificanceDF <- data.frame(ncol=2, nrow=length(unique(functionalSO\$term)))
count <- 1
for (functional in unique(functionalSO\$term)){
    sumTerms <- colSums(functionalSO[,3:6])
    functionalCount <- functionalSO[,3:6][functionalSO\$term==functional,]
    remainderCount <- sumTerms-functionalCount
    functionalMatrix <- as.matrix(rbind(functionalCount, remainderCount))
    # increase workspace to prevent fisher test from crashing
    # fisherFunctional <- fisher.test(functionalMatrix, workspace=100000000, hybrid=TRUE, simulate.p.value=TRUE)
    fisherFunctional <- chisq.test(functionalMatrix)
    functionalSignificanceDF[count,1] <- functional
    functionalSignificanceDF[count,2] <- fisherFunctional\$p.value
    count <- count +1
}


# Name columns
colnames(functionalSignificanceDF) <- c("terms","fisherPvalue")
# Perform Bonferonni adjustmet
functionalSignificanceDF\$fisherPvalueBonferroni <- p.adjust(functionalSignificanceDF\$fisherPvalue, method="bonferroni", n=length(functionalSignificanceDF\$fisherPvalue))
# See which classes differ significantly after multiple testing correction
functionalSignificanceDF[functionalSignificanceDF\$fisherPvalueBonferroni<0.05,]
functionalSignificanceDF[functionalSignificanceDF\$fisherPvalueBonferroni!=0,]

# Append p-values to frequency dataframe
functionalSOFreq <- merge(functionalSOFreq, functionalSignificanceDF, by="terms")

####### Order dataframe and make a index for plotting
# Add index for later plotting
functionalSignificanceDF\$index <- 1:nrow(functionalSignificanceDF)
significantLabels <- functionalSignificanceDF[functionalSignificanceDF\$fisherPvalueBonferroni<=0.05,]
significantLabels\$fisherPvalueBonferroni <- format(significantLabels\$fisherPvalueBonferroni,digits=3)

# Melt and reorder for plotting
mfunctionalSOFreq <- melt(functionalSOFreq, id=c("terms","fisherPvalueBonferroni","fisherPvalue" ))
# Remove combined allele frequency
mfunctionalSOFreq <- mfunctionalSOFreq[mfunctionalSOFreq\$variable!="freqAll",]
# Plot consequences that occur at a frequency of at least 1% in total
functionalFreqPlot <- ggplot(mfunctionalSOFreq, aes(x = terms, y = value) ) +
    geom_point(aes(shape=variable, col=variable), size=4) +
    theme_bw() + xlab("Consequence") +  ylab("Proportion") +
    scale_shape_manual(values=c(15:18)) +
    coord_cartesian(ylim = c(0.00, 0.45)) +
    scale_y_continuous(limits = c(0.00, 0.45)) +
    annotate("text", size=3, fontface="bold", y=0.19, angle=90,
    x=significantLabels\$index, label=significantLabels\$fisherPvalueBonferroni) +
    theme(legend.title=element_blank(), legend.position="bottom",legend.key=element_blank(),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face=c("bold")))
functionalFreqPlot
ggsave(dataset_pgxFunctionalClassesByFrequency_tiff, height=5.36*1.25, width=7.58, units='in', dpi=120)

# plot_grid(functionalSOplot, functionalFreqPlot,
# labels = c("A", "B"), align = "v", nrow=2)
# ggsave("pgxFunctionalClassesGrid.tiff",
# height=5.36*2, width=7.58, units='in', dpi=120)
