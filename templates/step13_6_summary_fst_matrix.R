#!/usr/bin/env Rscript

# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2018/03/04

library(ggplot2)

FstData <- read.table("AIBST.weighted.fst.estimates.matrix.txt", sep=',', h=T)
FstData <- as.data.frame(FstData)
FstMatrix <-data.matrix(FstData)

# Calculate summary of FST x 100
mean(FstMatrix,na.rm =T)
max(FstMatrix,na.rm =T)
min(FstMatrix,na.rm =T)