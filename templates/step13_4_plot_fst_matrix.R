#!/usr/bin/env Rscript

# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2018/03/04

library(ggplot2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

data <- read.table("${weighted_fst_estimates}")
names(data) <- c("POP1", "POP2", "FST")
upper_tri <- get_upper_tri(data)
# Create a ggheatmap

# Create a ggheatmap
ggheatmap <- ggplot(data, aes(POP2, POP1, fill = FST))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "red",
    midpoint = 0, space = "Lab",
    name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
    coord_fixed()

ggsave("${weighted_fst_estimates_png}", height=10, width=15.58, units='in', dpi=120)