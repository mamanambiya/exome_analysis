#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/08/02

## Create file with CHRM, POS, POP_1_FRQ, POP_2_FRQ

# Load Global allele frequency information
# myFiles <- list.files(path = myPath, pattern="*.Ind.frq")

## These variables are from nextflow
GROUP_POP = "${dataset}"
myFiles = as.list(strsplit("${POP_Ind_frq_files}", split=" "))[[1]]
populations = as.list(strsplit("${POPs}", " "))[[1]]

for (file in myFiles){
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
        # Still need to remove first population data
        dataset <- read.table(file, header = T, sep = "\\t")
        dataset <- dataset[, 1 : 3]
    }
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
        temp_dataset <- read.table(file, header = T, sep = "\\t")
        dataset <- cbind(dataset, temp_dataset[, 4])
        colnames(dataset)[length(dataset)] <- colnames(temp_dataset)[length(temp_dataset)]
        rm(temp_dataset)
    }
}

## Add column names
# popColnames <-  c("CHROM", "POS")
# popColnames <-  c("rsID", "REF", "ALT")
# for (column in (length(dataset2)-3)) {
#     for (i in 1:length(populations)){
#         popFreq <- paste0(populations[i], "_Alt_FREQ")
#         popColnames <- append(popColnames, popFreq)
#     }
# }
# colnames(dataset) <- popColnames

# Upload info and rsID file
INFOandrsID <- read.table("${GROUP_biall_info_rsid}", header = T, sep = "\\t")
# Combine annotated with dataset, removing duplicated first two columns of dataset
# datasetAnnotated <- cbind(INFOandrsID, dataset[,-(2:3)])
datasetAnnotated <- cbind(INFOandrsID[, - (12 : 18)], dataset[, - (2 : 3)], INFOandrsID[, (12 : 18)])

# # If mean allele frequency is > 0.5 flip alleles in those columns
# for (row in 1:nrow(datasetAnnotated)) {
#     if (datasetAnnotated\$AF[row]>0.5) {
#         datasetAnnotated[row,grep('*FREQ',names(datasetAnnotated))] <-
#         lapply(datasetAnnotated[row,grep('*FREQ',names(datasetAnnotated))], function(x) 1-x)
#     }
# }

### GET NUMBER of POLYMORPHIC ALLELES IN EACH POPULATION
datasetFreq <- datasetAnnotated[,grep('*FREQ',names(datasetAnnotated))]
polymorphic <- colSums(datasetFreq != 0)
write.csv(polymorphic, paste(GROUP_POP, "_dataset_polymorphic.csv", sep=''), row.names = FALSE)


### GET NUMBER of SINGLETON ALLELES IN EACH POPULATION

# Create unique ID column for matching
datasetAnnotated\$UniqueID <- paste(datasetAnnotated\$CHROM, datasetAnnotated\$POS, datasetAnnotated\$ALT, sep=":")

# Add a column indicating whether the variant is a singleton or not
datasetAnnotated\$Singleton <- with(datasetAnnotated, ifelse(AC==1, "Yes","No"))
# Add a column indicating whether the variant has a global MAF>=0.01
datasetAnnotated\$Maf0.01 <- with(datasetAnnotated, ifelse(AF>=0.1, "Yes","No"))

## Save dataset and datasetAnnotated into files for future use
write.csv(dataset, paste(GROUP_POP, "_dataset.csv", sep=''), row.names = FALSE)
write.csv(datasetAnnotated, paste(GROUP_POP, "_datasetAnnotated.csv", sep=''), row.names = FALSE)
