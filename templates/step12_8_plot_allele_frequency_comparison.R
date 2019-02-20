#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(ggplot2)
library(reshape2)
library(plyr)

###
### Step 12.8: Plot Allele frequency comparison
###

# Plot and calculate R2
# dataset1 <- "${dataset1}"
# dataset2 <- "${dataset2}"
datasets_data <- read.table("${datasets_frq}", header=T)
colnames(datasets_data) <- c("CHRM:POS", "D2", "D1")
datasets_R2 <- round((summary(lm(datasets_data\$D1 ~ datasets_data\$D2))\$r.squared),digits = 2)
datasets_Plot <-ggplot(data=datasets_data, aes(x=D1, y=D2)) + geom_point() +
  theme_bw() + xlab("${dataset1} frequency") + ylab("${dataset2} frequency") +
  geom_smooth(method='lm') +
  annotate("text", size=5, x = 0.1, y = 0.95, label = paste("R2 =",datasets_R2))

ggsave("${ExternalAlleleFreqComparison_tiff}",
       height=5.36, width=7.58, units='in', dpi=120)
