library(Seurat)
library(rhdf5)
library(RANN)
library(tidyr)
library(dplyr)
library(SingleCellExperiment)
library(Matrix)
library(limma)
library(edgeR)

source("scratch/Group4/github_original/first_pipeline/our_functions.R")

#cells <- c(100, 200, 600, 1000, 1500, 2000, 4000, 6000, 8000)
cells <- c(100, 2000, 4000, 6000, 8000)
umis <- c(10000, 2000, 300, 50, 10)
d <- "bipolar"
for (c in rev(cells)){
  for (u in umis){
    print(c)
    print(u)
    loom_file <- sprintf("/home/jovyan/scratch/Group4/subsampled/%s_cells_%s_umis_%s.loom", d, c, u)
    sce <- loom2sce(loom_file)
    sce <- processData(sce)
    clusters <- levels(colData(sce)$cluster_ids)
    all <- lapply(1:(length(clusters)-1), function(i){
      sub <- lapply((i+1):length(clusters), function(j){
        ourFeatures(sce, clusters[i], clusters[j])
      })
      do.call(rbind, sub)
    })
    allFeat <- do.call(rbind, all)
    allFeat <- as.data.frame(allFeat)
    allFeat[, 3:ncol(allFeat)] <- apply(allFeat[, 3:ncol(allFeat)], 2, as.numeric)
    write.csv(allFeat, paste0(gsub('.loom', '_features.csv', loom_file)))
  }
}


loom_file <- "/home/jovyan/data/tasks/how_many_cells/bipolar.loom"
sce <- loom2sce(loom_file)
sce <- processData(sce)
clusters <- levels(colData(sce)$cluster_ids)
logFClist <- lapply(1:(length(clusters)-1), function(i){
  sub <- lapply((i+1):length(clusters), function(j){
    # cell names for clusters
    c1type <- clusters[i]
    c2type <- clusters[j]
    c1names <- rownames(colData(sce))[colData(sce)$cluster_ids == c1type]
    c2names <- rownames(colData(sce))[colData(sce)$cluster_ids == c2type]
    
    # means
    s1 <- Matrix::rowMeans(assay(sce)[, c1names])
    s2 <- Matrix::rowMeans(assay(sce)[, c2names])
    
    #logFC
    dge <- DGEList(counts=data.frame(c1 = s1, c2 = s2))
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    
    c(c1 = c1type, c2 = c2type, logFC = var(logCPM[, 1] - logCPM[, 2]))
  })
  do.call(rbind, sub)
})
logFC <- as.data.frame(do.call(rbind, logFClist))
logFC[, 3] <- as.numeric(as.vector(logFC[, 3]))
write.csv(logFC, "/home/jovyan/scratch/Group4/subsampled/bipolar_logFC.csv")


