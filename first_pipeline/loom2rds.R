#Rscript loom2rds.R
#this converts the loom file into single cell experiment object, including PCA factors
pth<-"/home/jovyan/scratch/Group4"
spth<-file.path(pth,"github_original/first_pipeline")
ipth<-file.path(pth,"subsampled")
opth<-file.path(pth,"subsampled_sce")
source(file.path(spth,"our_functions.R"))
dir.create(opth)
#loom_file <- "/home/jovyan/data/tasks/how_many_cells/human_pancreas.loom"
#iname<-"pbmc_cells_2000_umis_300.loom"
f<-function(iname){
  print(paste("Now processing:",iname))
  sce<-loom2sce(file.path(ipth,iname))
  sce<-processData(sce)
  oname<-strsplit(iname,".",fixed=TRUE)[[1]][1]
  saveRDS(sce,file=file.path(opth,paste0(oname,".rds")))
}
res<-mclapply(list.files(ipth),f,mc.cores=15)

