#!/usr/bin/env Rscript
#  Title     : Plot population frequencies supergroup population
# Objective : Rare Global Common Variants
# Created by: mamana
# Created on: 2017/08/11

library(ggplot2)
library(reshape2)

###
### Step 8
####

superpop.plot.colours <- c("#009E73","#CC79A7", "#56B4E9","#E69F00","#000000")

# ADD POPULATION DATA to the DATAFRAME
samples <- read.table("${GROUP_POP_sample}", header = F)
colnames(samples) <- c("sample", "pop", "supergroup","sex", "dataset")
# Make data frame with unique populations and their superpops
popDetail <- unique(samples[,2:3])
colnames(popDetail) <- c("pop", "supergroup")
popDetail <- popDetail[with(popDetail, order(popDetail\$pop)), ]
populations <- popDetail\$pop

# Load annotated data
datasetAnnotated <- read.table("${GROUP_POP_datasetAnnotated}", sep=',', header=T)
# Add CADD scores
cadd <- read.table("${cadd_annotations}")
colnames(cadd) <- c("CHROM", "POS", "REF", "ALT", "RawScore", "PHRED")
cadd\$UniqueID <- paste(cadd\$CHROM, cadd\$POS, cadd\$ALT, sep=":")
cadd <- merge(cadd, datasetAnnotated, by="UniqueID")
cadd.delterious <- cadd[cadd\$PHRED>=20,]

# Plot highly differentiated SNPs (5% in one population, but 0.5% in the global population)
highDiffPop <- read.table("${GROUP_POP_population_pop_diff_snps_genotypes}", header=F)
colnames(highDiffPop) <- c("Population", "ID")
highDiffGenes <- read.table("${GROUP_POP_genes_pop_diff_snps_genotypes}", header=F)
colnames(highDiffGenes) <- c("Gene", "Function", "ID", "UniqueID")
highDiffGenes\$Function <- gsub("_", " ", highDiffGenes\$Function)

# Write table of population differentiated SNPs
## highDiff <- merge(highDiffPop, highDiffGenes, by="ID") ## Not working because of size. use Data.table library

highDiff <- read.table("${GROUP_POP_population_gene_pop_diff_snps_genotypes}", header=F)
colnames(highDiff) <- c("Population", "Gene", "Function", "ID", "UniqueID")
highDiff <- merge(highDiff, cadd[,c(1,7)], by = "UniqueID")
highDiff <- highDiff[order(highDiff\$Gene),]
write.csv(highDiff, "${GROUP_POP_csv}", row.names = FALSE, quote = FALSE)

# Make dataset for plotting
highDiffrsIDFreq <- merge(highDiff, datasetAnnotated, by="UniqueID")
highDiffrsIDFreq\$Gene <- paste("(", highDiffrsIDFreq\$Gene, ":", sep ="")
highDiffrsIDFreq\$ID <- paste(highDiffrsIDFreq\$ID, ")", sep ="")
highDiffrsIDFreq\$Label <- paste(highDiffrsIDFreq\$Population, highDiffrsIDFreq\$Gene, highDiffrsIDFreq\$ID, sep = " ")
highDiffrsIDFreq <- highDiffrsIDFreq[, c(grep('*Label*', names(highDiffrsIDFreq)), grep('*_Alt_FREQ', names(highDiffrsIDFreq)))]

# Make tall data frame for plotting in ggplot
mHighDiffrsIDFreq <- melt(highDiffrsIDFreq, id="Label")
mHighDiffrsIDFreq\$Pop <- rep(sort(popDetail\$pop), each=nrow(highDiffrsIDFreq))
mHighDiffrsIDFreq\$Super <- rep(popDetail\$supergroup[order(popDetail\$pop)],each=nrow(highDiffrsIDFreq))

# Plot allele frequency stripchart
# Allele frequencies of highly differentiated pharmacogenomic variants
highDiffPlot <- ggplot(mHighDiffrsIDFreq, aes(x = Label, y = value) ) +
    geom_point( aes(colour = pop), shape=19) +
    theme_bw() + xlab("Variant") +
    scale_colour_manual(values=superpop.plot.colours) +
    ylab("Frequency") + xlab("Population and variant") + geom_hline(aes(yintercept=0.05), linetype="dashed") +
    theme(legend.title=element_blank(),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="bold"))
highDiffPlot
ggsave("${GROUP_POP_tiff}", height=5.36, width=7.58, units='in', dpi=120)
