#Rscript /home/jovyan/scratch/Group4/github_original/first_pipeline/diff_expr.R
library(SingleCellExperiment)
library(edgeR)
library(limma)

diff_expr<-function(cids,m,cm,alpha=.05){
  #alpha=threshold for FDR significance
  #m a count matrix
  #cids a set of cluster ids to subset by
  #cm colData corresponding to count matrix
  #returns:
  #1. the variance of the log fold changes across all genes
  #2. Number of "differentially expressed" genes at significance (FDR) alpha
  grp<-factor(cm$cluster_ids[cm$cluster_ids %in% cids])
  m<-m[,cm$cluster_ids %in% cids]
  gg<-Matrix::rowSums(m)>0 #remove genes that are all zero
  #naive application of limma-trend
  dge <- DGEList(counts=m[gg,])
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  X<-model.matrix(~grp)
  fit <- lmFit(logCPM, design=X)
  fit <- tryCatch(eBayes(fit, trend=TRUE),error=function(e){NULL})
  if(is.null(fit)){
    return(c(NA,NA))
  } else {
    tt<-topTable(fit,coef=ncol(X),number=nrow(logCPM))
    #lfcvar<-var(tt$logFC)
    #tt<-tt[tt$adj.P.Val<alpha,]
    #p05<-nrow(res)
    #minLFC<-min(abs(res$logFC))
    #with(res,plot(logFC,-log10(P.Value))) #volcano plot
    return(c(var(tt$logFC),sum(tt$adj.P.Val<alpha)))
  }
}

diff_expr_all<-function(iname,pth="/home/jovyan/scratch/Group4",verbose=FALSE){
  print(paste("Now processing:",iname))
  ipth<-file.path(pth,"subsampled_sce")
  opth<-file.path(pth,"subsampled_dge")
  dir.create(opth)
  #iname<-"bipolar_cells_2000_umis_2000.rds"
  #iname<-"pbmc_cells_8000_umis_2000"
  sce<-readRDS(file.path(ipth,iname))
  cm<-colData(sce)
  cid_uniq<-levels(cm$cluster_ids)
  m<-assay(sce,"counts")
  res<-data.frame(matrix(NA,nrow=choose(length(cid_uniq),2),ncol=4))
  colnames(res)<-c("cl1","cl2","var_log_fc","num_de_genes")
  t<-0
  for(i in 1:(length(cid_uniq)-1)){
    for(j in (i+1):length(cid_uniq)){
      t<-t+1
      if(verbose) print(paste(t,"/",nrow(res)))
      cids<-cid_uniq[c(i,j)]
      res[t,1:2]<-cids
      res[t,3:4]<-diff_expr(cids,m,cm)
    }
  }
  oname<-strsplit(iname,".",fixed=TRUE)[[1]][1]
  res$file<-oname
  write.csv(res,file=file.path(opth,paste0(oname,".csv")),quote=FALSE,row.names=FALSE)
}

ipth<-"/home/jovyan/scratch/Group4/subsampled_sce"
q<-mclapply(list.files(ipth),diff_expr_all,verbose=FALSE,mc.cores=15)



