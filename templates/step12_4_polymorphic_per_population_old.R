#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(ggplot2)
library(reshape2)
library(plyr)

###
### Step 12.4: Plot singleton per population
###

## ANALYSE NUMBER of SINGLETONS/POPULATION and PLOT TOTAL POLYMORPHIC/POPULATION

# wdir <- "/home/mamana/Project_data/exome_aibst/data/MERGED_POP/FRQ/AIBST"
# ADD POPULATION DATA to the DATAFRAME
# setwd(wdir)
# Make superpop palette for x axis labels
superpopXcolours <- c(rep("#009E73",7), rep("#CC79A7",4), rep("#56B4E9",5), rep("#E69F00",5), rep("#000000",5))
# Make two colour (colour blind friendly) palette for two colour plots
two.colour.palette <- c("#CC79A7", "#0072B2")
samples <- read.table("${dataset_sample}",header = F)
colnames(samples) <- c("sample", "pop", "supergroup","sex", "dataset")
# Make data frame with unique populations and their superpops
popDetail <- unique(samples[,2:3])
colnames(popDetail) <- c("pop", "supergroup")
popDetail <- popDetail[with(popDetail, order(popDetail\$pop)), ]
populations <- popDetail\$pop
# populations <- append(as.vector(popDetail\$pop), "KG_AFR")
# populations <- append(populations, "gnomAD_AFR")
# populations <- append(populations, "AGVP")

datasetAnnotated <- read.table("${dataset_datasetAnnotated}", header=T, sep=',')
# ### GET NUMBER of POLYMORPHIC ALLELES IN EACH POPULATION
datasetFreq <- datasetAnnotated[,grep('*FREQ',names(datasetAnnotated))]
polymorphic <- colSums(datasetFreq != 0)
#polymorphic
datasetFreqSing <- datasetAnnotated[datasetAnnotated\$Singleton=="Yes",c(grep('*_Alt_FREQ',names(datasetAnnotated)))]
# Count number of non-zero occurences in each
popSingles <- apply(datasetFreqSing, 2, function(x) sum(x!=0))
polymorphicTable <- cbind(polymorphic, popSingles)
rownames(polymorphicTable) <- populations
colnames(polymorphicTable) <- c("Polymorphic", "Singletons")
polymorphicTable <- as.data.frame(polymorphicTable)
polymorphicTable\$NonSingleton <- polymorphicTable\$Polymorphic - polymorphicTable\$Singletons
polymorphicTable\$Populations <- populations
polymorphicTable\$Polymorphic <- NULL
polymorphicTableM <- melt(polymorphicTable, id.vars = "Populations")
# Reorder for plotting
polymorphicTableM\$Populations <- factor(polymorphicTableM\$Populations,levels=popDetail\$pop[order(popDetail\$pop)])

polymorphicTableM\$variable <- gsub("NonSingleton", "Non-singletons", polymorphicTableM\$variable)
# Plot singletons vs non-singletons per population using a stacked bar graph
ggplot(polymorphicTableM, aes(x=Populations, y=value, fill=variable))  +
  geom_bar(stat="identity") + xlab("Population") + ylab("Count") + theme_bw() +
  theme(legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=12,face="bold",
   colour=superpopXcolours)) + scale_fill_manual(values=two.colour.palette)
ggsave("${dataset_pgxPolymorphicPerPopulation}", height=10, width=15.58, units='in', dpi=120)
