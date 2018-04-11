library(plyr)
library(stringr)
#merge data
pth<-"/home/jovyan/scratch/Group4"
f<-function(x,pth){
  read.csv(file.path(pth,x),header=TRUE)
}
dpth<-file.path(pth,"subsampled_dge")
dge_results<-do.call("rbind",lapply(list.files(dpth),f,pth=dpth))
fnames<-as.character(unique(dge_results$file))
spth<-file.path(pth,"subsampled")
g<-function(x,pth,verbose=FALSE){
  if(verbose) print(x)
  csv_file<-file.path(pth,paste0(x,"_features.csv"))
  if(file.exists(csv_file)){
    q<-read.csv(csv_file)
    q<-q[,-1]
    if(is.null(ncol(q))) return(NULL)
    q$file<-x
    colnames(q)[1:2]<-c("cl1","cl2")
    return(q)
  } else {
    if(verbose) print(paste(csv_file,"not found"))
    return(NULL)
  }
}
sep_results<-do.call("rbind",lapply(fnames,g,pth=spth,verbose=FALSE))
#x<-fnames[1]
#q<-read.csv(file.path(spth,paste0(x,"_features.csv")),header=TRUE)
res<-join(sep_results,dge_results,c("cl1","cl2","file"),type="inner")
colnames(res)[colnames(res)=="ncells2"]<-"ncell2"
fname2cells_umi<-function(fnames){
  #q<-matrix(unlist(strsplit(fnames,"_",fixed=TRUE)),nrow=5)
  q<-matrix(as.numeric(unlist(str_extract_all(fnames,"\\d+"))),nrow=2)
  data.frame(cells=q[1,],umis=q[2,])
}
#cu<-fname2cells_umi(res$file)
res<-cbind(res,fname2cells_umi(res$file))
res$c1freq<-with(res,ncell1/cells)
res$c2freq<-with(res,ncell2/cells)
write.csv(res,file=file.path(pth,"merged_features.csv"),row.names=FALSE,quote=FALSE)
