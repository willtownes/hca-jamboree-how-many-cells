loom2sce<-function(loom_file){
  #read a loom_file with single cell expression data from file path loom_file
  #returns a SingleCellExperiment object with sparse Matrix expression values
  #coldata includes cluster labels and sequencing depth for each cell
  data <- h5read(loom_file,name = "/")
  data_matrix <- t(data$matrix)
  rownames(data_matrix)<-as.character(data$row_attrs$gene_names)
  colnames(data_matrix)<-as.character(data$col_attrs$cell_names)
  cluster_ids <- as.character(data$col_attrs$cluster)
  rm(data)
  data_matrix <- Matrix(data_matrix,sparse=T)
  cmeta<-data.frame(cluster_ids=cluster_ids,umi_depth=Matrix::colSums(data_matrix))
  stopifnot(all(colnames(data_matrix)==rownames(cmeta)))
  res<-SingleCellExperiment(assays=list(counts=data_matrix))
  colData(res)<-DataFrame(cmeta)
  res
}


processData <- function(sce, npca = 30) {
  seurat <- CreateSeuratObject(raw.data = assay(sce))
  normSeurat <- NormalizeData(seurat,
                              display.progress = F,
                              scale.factor = 1e4,
                              normalization.method = "LogNormalize")
  varGenes <- FindVariableGenes(normSeurat,
                                do.plot = FALSE,
                                display.progress = FALSE)
  varGenes@var.genes <- rownames(head(varGenes@hvg.info, 1000))
  scale <- ScaleData(varGenes, 
                     genes.use = varGenes@var.genes,
                     display.progress = FALSE,
                     check.for.norm = FALSE)
  pca <- RunPCA(scale,
                pcs.compute = npca,
                do.print = FALSE)
  pca <- GetCellEmbeddings(pca, reduction.type = "pca", dims.use = 1:npca)
  reducedDims(sce) <- SimpleList(PCA=pca)
  sce
}

computeSeparability <- function(input.data, cells.1, cells.2, k = 20) {
  if (length(cells.1) < 3) {
    #stop("Fewer than 3 cells in the first group (cells.1)")
    return(NA)
  }
  if (length(cells.2) < 3) {
    #stop("Fewer than 3 cells in the second group (cells.2)")
    return(NA)
  }
  k <- min(c(k, length(cells.1), length(cells.2)))
  tnn <- nn2(data = input.data[c(cells.1, cells.2), ], k = k + 1)
  idx <- tnn$nn.idx[, -1]
  rownames(idx) <- c(cells.1, cells.2)
  
  N1 <- length(cells.1)
  N2 <- length(cells.2)
  correct_neighbors_c1 <- rowSums(idx[seq_len(N1),] <= N1)
  correct_neighbors_c2 <- rowSums(idx[seq_len(N2) + N1,] > N1)
  
  return(mean(c(correct_neighbors_c1, correct_neighbors_c2)) / k)
}

runLimmavoom <- function(sce) {
  colData(sce)[, 'cluster_ids'] <- droplevels(colData(sce)[, 'cluster_ids'])
  design <- model.matrix(~ colData(sce)[, 'cluster_ids'])
  dgel <- DGEList(assay(sce))
  v <- voom(dgel, design, plot=FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, n = nrow(dgel), sort.by = "none")
  pvals <- tt$P.Value
  padj <- p.adjust(pvals, method = "BH")
  padj[is.na(padj)] <- 1
  data.frame(gene = rownames(tt), pval=pvals, padj=padj, logfc=tt$logFC)
}

compute_n_de <- function(sce, alpha = .05){
  cl <- factor(colData(sce)$cluster_ids)
  gg <- Matrix::rowSums(assay(sce)) > 0 #remove genes that are all zero
  #naive application of limma-trend
  dge <- DGEList(counts = assay(sce)[gg, ])
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log = TRUE, prior.count = 3)
  X <- model.matrix(~ cl)
  fit <- lmFit(logCPM, design = X)
  fit <- tryCatch(eBayes(fit, trend = TRUE), error = function(e){NULL})
  if (is.null(fit)){
    nDE <- NA
  } else {
    tt <- topTable(fit, coef = ncol(X), number = nrow(logCPM))
    nDE <- sum(tt$adj.P.Val < alpha)
  }
  return(nDE)
}

compute_var_logFC <- function(sce, c1names, c2names){
  s1 <- Matrix::rowMeans(assay(sce)[, c1names])
  s2 <- Matrix::rowMeans(assay(sce)[, c2names])
  dge <- DGEList(counts=data.frame(c1 = s1, c2 = s2))
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  varlogFC <- var(logCPM[, 1] - logCPM[, 2])
  return(varlogFC)
}


ourFeatures <- function(sce, c1type, c2type, k = 20, thr = 0.05){
  # cell names for cell type 1 and 2
  c1names <- rownames(colData(sce))[colData(sce)$cluster_ids == c1type]
  c2names <- rownames(colData(sce))[colData(sce)$cluster_ids == c2type]
  
  # separability
  cluster_sep <- computeSeparability(reducedDim(sce), c1names, c2names, k = k)
  
  meanLibSize = mean(Matrix::colSums(assay(sce)))
  if (class(assay(sce)[, c1names]) == "dgCMatrix"){
    c1LibSize = mean(Matrix::colSums(assay(sce)[, c1names]))
  }else{
    c1LibSize = sum(assay(sce)[, c1names]) 
  }
  if (class(assay(sce)[, c2names]) == "dgCMatrix"){
    c2LibSize = mean(Matrix::colSums(assay(sce)[, c2names]))  
  }else{
    c2LibSize = sum(assay(sce)[, c2names])  
  }
  
  # return
  c(c1 = c1type,
    c2 = c2type,
    ncell1 = length(c1names),
    ncell2 = length(c2names),
    sep = cluster_sep,
    meanLibSize = meanLibSize,
    c1LibSize = c1LibSize,
    c2LibSize = c2LibSize
    )
}
