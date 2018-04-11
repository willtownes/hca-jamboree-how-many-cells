#!/usr/bin/env Rscript

suppressMessages({
  library(Seurat)
  library(rhdf5)
  library(RANN)
  library(tidyr)
  library(dplyr)
  library(SingleCellExperiment)
  library(Matrix)
  library(limma)
  library(edgeR)
})

source("~/scratch/Group4/hca-jamboree-how-many-cells/first_pipeline/our_functions.R")

# from python train model
paramEstimates <- c(intercept = 0.7072,
                    cfreq = 0.4125,
                    ncell = 0.5842,
                    meanLibSize = -.1163,
                    varLogFC = 0.3915,
                    nDE = 0.3079)

inv_logit <- function(x){
  1/(1 + exp(-x))
}

make_compute_separability <- function(paramEstimates){
  fct <- function(ncell, meanLibSize, fcluster1, fcluster2, varLogFC, nDE){
    x <- sum(paramEstimates * c(1, log(fcluster1 + fcluster2), 
                                log(ncell), log(meanLibSize), 
                                log(varLogFC), log(nDE)))
    pred <- inv_logit(x)
    return(pred)
  }
  return(fct)
}

# read args
args <- commandArgs(trailingOnly = TRUE)
loom_file <- args[1]
fcluster1 <- as.numeric(args[2])
fcluster2 <- as.numeric(args[3])
tsvFile <- args[4]

run = FALSE 
if (run){
  loom_file <- 'scratch/Group4/subsampled/human_pancreas_cells_2000_umis_2000.loom'
  fcluster1 <- 0.03
  fcluster2 <- 0.02
  tsvFile <- 'mypred.csv'
}

# loom format to singleCellExperiment object
print(sprintf("Reading loom file %s", loom_file))
print("Converting into SingleCellExperiment")
sce <- loom2sce(loom_file)

# normalize and perform pca 
print("Pre-processing")
sce <- processData(sce)

# mean lib size
print("Compute mean library size")
## overall 
meanLibSize <- mean(Matrix::colSums(assay(sce)))

## cluster1
print("Compute mean library size for cluster1")
clusters <- levels(factor(colData(sce)$cluster_ids))
c1type <- clusters[1]
c1names <- rownames(colData(sce))[colData(sce)$cluster_ids == c1type]
c1LibSize <- mean(Matrix::colSums(assay(sce)[, c1names]))

## cluster2
print("Compute mean library size for cluster2")
c2type <- clusters[2]
c2names <- rownames(colData(sce))[colData(sce)$cluster_ids == c2type]
c2LibSize <- mean(Matrix::colSums(assay(sce)[, c2names]))  

# variance of log fold change
print("Compute variance of logFC")
varLogFC <- compute_var_logFC(sce, c1names, c2names)

# nDE
print("Compute nDE")
nDE <- compute_n_de(sce)


# prediction of separability as a function of ncell
ncell <- seq(1000, 100000, 1000)
print("Predict separability")
compute_separability <- make_compute_separability(paramEstimates)
predSep <- sapply(ncell, 
                  compute_separability,
                  meanLibSize, fcluster1, fcluster2, varLogFC, nDE)


# write predicted separability in tsv file 
print(sprintf("Write predicted separability to %s", tsvFile))
predDf <- data.frame(ncell = ncell, predSeparability = predSep)
colnames(predDf) <- c("Number of Cells", "Predicted Separability")
write.table(x = predDf, file = tsvFile, sep = "\t", quote = FALSE, row.names = FALSE)



