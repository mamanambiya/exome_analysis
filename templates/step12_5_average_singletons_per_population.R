#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(ggplot2)
library(reshape2)
library(plyr)

## Step 12.5

dataset_sample <- "${dataset_sample}"
dataset_singletons_per_sample <- "${dataset_singletons_per_sample}"
dataset_pgxAverageSingletonsPerPopulation <- "${dataset_pgxAverageSingletonsPerPopulation}"

# Make superpop palette for x axis labels
superpopXcolours <- c(rep("#009E73",7), rep("#CC79A7",4), rep("#56B4E9",5), rep("#E69F00",5), rep("#000000",5))
# Make two colour (colour blind friendly) palette for two colour plots
two.colour.palette <- c("#CC79A7", "#0072B2")
samples <- read.table(dataset_sample,header = F)
colnames(samples) <- c("sample", "pop", "supergroup","sex", "dataset")
# Make data frame with unique populations and their superpops
popDetail <- unique(samples[,2:3])
colnames(popDetail) <- c("pop", "supergroup")
popDetail <- popDetail[with(popDetail, order(popDetail\$pop)), ]
populations <- popDetail\$pop

# Calculate number of singletons / individual / population
popSinglesCount <- read.table(dataset_singletons_per_sample, header =T)
popSinglesCount <- popSinglesCount[,-c(3,4)]

popSinglesSummary <- ddply(popSinglesCount, "pop", summarise,
    mean = mean(singletons), sd = sd(singletons),
    sem = sd(singletons)/sqrt(length(singletons)))

popSinglesSummary\$pop <- factor(popSinglesSummary\$pop,levels=popDetail\$pop[order(popDetail\$supergroup)])

ggplot(popSinglesSummary, aes(x=pop, y=mean))  +
    geom_bar(stat="identity", fill="#009E73") + xlab("Population") +
    ylab("Count") + theme_bw() +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem)) +
    theme(legend.title=element_blank(),
    axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=12,
    face="bold", colour=superpopXcolours))
ggsave(dataset_pgxAverageSingletonsPerPopulation, height=10, width=15.58, units='in', dpi=120)