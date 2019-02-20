#!/usr/bin/env Rscript

# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2019/01/18


library(plyr)
library(ggplot2)
library(reshape2)

superpop.plot.colours <- c("#3A8B8B","#7CFC00","#CD9B1D","#7FFFD4","#FF8C00","#6495ED","#008B8B","#FF0000", "#008B45","#8B4500","#8B7355","#CDCD00","#98FB98","#8B4513","#FFD700","#00FFFF", "#98F5FF","#228B22","#548B54","#8B3A62","#CD6090","#CDC673","#FF69B4","#FFEC8B", "#CD6600","#2E8B57")

# Analyse highly differentiated polymorphisms
pgxVariants <- read.table("${fst_csv}", header=T, fill=T, sep=',')

# Replace NAs with 0
pgxVariants[is.na(pgxVariants)] <- 0

# Select columns where at least one of the Fst comparisons
# is highly differentiated (i.e. >0.5)

pgxVariantsHighFst_ann <- read.table("${fst_HighDiff_ann}", header=T, fill=T, sep='\\t')
pgxVariantsHighFst_ann[is.na(pgxVariantsHighFst_ann)] <- 0

pgxVariantsHighFst_nann <- read.table("${fst_HighDiff_csv}", header=T, fill=T, sep=',')
pgxVariantsHighFst_nann[is.na(pgxVariantsHighFst_nann)] <- 0
pgxVariantsHighFst_nann <- rename(pgxVariantsHighFst_nann, c("CHRM.POS"="CHROM.POS"))

# Calculate average Fst for each variant
pgxVariantsHighFst_nann\$Fst.mean <- rowMeans(pgxVariantsHighFst_nann[,grep('*Fst',names(pgxVariantsHighFst_nann))])
pgxVariantsHighFst_nann\$Fst.mean <- pgxVariantsHighFst_nann\$Fst.mean
row.names(pgxVariantsHighFst_nann) <- NULL

# Add Ethnic groups freq info
datasetAnnotated <- read.table("${datasetAnnotated}", header = T, sep = ",")
datasetAnnotated\$CHROM.POS <- paste(datasetAnnotated\$CHROM, datasetAnnotated\$POS, sep=":")

# Merge pgxVariantsHighFst_ann and pgxVariantsHighFst_nann
pgxVariantsHighFst_ <- merge(pgxVariantsHighFst_ann, pgxVariantsHighFst_nann, by="CHROM.POS")
pgxVariantsHighFst <- merge(pgxVariantsHighFst_, datasetAnnotated, by="CHROM.POS")

pgxVariantsHighFst <- rename(pgxVariantsHighFst, c("ANN.0..EFFECT"="FUNCTION", "ANN.0..GENE"="GENE"))
# Remove function aliases from those that have them
pgxVariantsHighFst\$FUNCTION <- gsub("&.*","",pgxVariantsHighFst\$FUNCTION)

# Select rows with the maximum mean FST per gene
pgxVariantsHighFstMax <- ddply(pgxVariantsHighFst, .(GENE), function(x) x[which.max(x\$Fst.mean),])
# Remove multiallelic
# pgxVariantsHighFstMax <- pgxVariantsHighFstMax[-grep(",", pgxVariantsHighFstMax\$ALT),]
# Add CADD scores
#pgxVariantsHighFstMax <- merge(pgxVariantsHighFstMax, cadd[,c(7,15)])
#pgxVariantsHighFstMax <- pgxVariantsHighFstMax[,-c(2:7)]
pgxVariantsHighFstMax <- pgxVariantsHighFstMax[order(pgxVariantsHighFstMax\$GENE),]
write.csv(pgxVariantsHighFstMax, "${HighGene_out}", row.names = FALSE)

# Add frequencies for these polymorphisms for plotting
#highFstDataset <- merge(pgxVariantsHighFstMax, datasetAnnotated, by="ID")
highFstDataset <- pgxVariantsHighFstMax

# Select all columns with ref allele frequency
RefhighFstDataset <- highFstDataset[,c(grep('*_Alt_FREQ',names(highFstDataset)))]
#RefhighFstDataset <- highFstDataset
RefhighFstDataset\$ID <-  sapply(strsplit(as.character(highFstDataset\$ID.x), ";"), '[',1)
RefhighFstDataset\$Gene <- highFstDataset\$GENE
RefhighFstDataset\$Function <- as.factor(highFstDataset\$FUNCTION)
RefhighFstDataset\$Fst.mean <- highFstDataset\$Fst.mean
RefhighFstDataset\$GeneRs <- paste(RefhighFstDataset\$Gene, RefhighFstDataset\$ID, sep=" ")
# RefhighFstDataset\$GeneRs <- paste(RefhighFstDataset\$Gene)
RefhighFstDataset <- RefhighFstDataset[order(RefhighFstDataset\$Gene),]

# Plot allele frequencies
# Make tall data frame for plotting in ggplot
mRefhighFstDataset <- melt(RefhighFstDataset, id=c("ID","Gene","Function","Fst.mean","GeneRs"))
mRefhighFstDataset\$variable <- gsub("_Alt_FREQ", "", mRefhighFstDataset\$variable)
# mRefhighFstDataset\$Pop <- rep(sort(popDetail\$pop),each=nrow(RefhighFstDataset))
# mRefhighFstDataset\$Super <- rep(popDetail\$super_pop[order(popDetail\$pop)],each=nrow(highFstDataset))

# Reorder factors for later plotting in ggplot from max mean FST
mRefhighFstDataset\$GeneRs <- factor(mRefhighFstDataset\$GeneRs,levels=RefhighFstDataset[order(RefhighFstDataset\$Fst.mean),"GeneRs"])

# Plot allele frequency stripchart
# Allele frequencies of highly differentiated pharmacogenomic variants
ggplot(mRefhighFstDataset, aes(x = GeneRs, y = as.numeric(value)) ) +
    geom_point(aes(colour=variable), shape=19, size=1) +
    theme_bw() + xlab("Variant") +
    scale_colour_manual(values=superpop.plot.colours) +
    ylab("Frequency") +
    scale_y_continuous() +
    theme(legend.title=element_blank(),
    axis.text.x=element_text(angle = 90, hjust = 1,
    vjust = 0.5, size=10, face="bold"))

ggsave("${HighFst_plot}",
height=6, width=25, units='in', dpi=120)

# scale_shape_manual(values=ifelse(variable == "AIBST", c(15,16,17) ,c(0,1,2)) , guide = "none")
# scale_shape_manual(values = c('AIBST' = 17, 'Men' = 16))