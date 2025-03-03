#!$HOME/bin/Rscript

HOME <- getwd()

args <- commandArgs(T)

lib <- args[1]
type <- args[2]

cov.list <- list()

chr <- c("chrI","chrII","chrIII","chrIV","chrV","chrX")

for (x in chr){
  cov.list[[x]] <- read.table(file=paste(HOME,"/",lib,"/",lib,".",x,"_",type,"_rpm.txt",sep=""),header=F,sep="\t")
}

cov.df <- data.frame(NULL)
for (y in 1:length(chr)){
  cov.df <- rbind.data.frame(cov.df,cov.list[[y]])
}

write.table(x=cov.df,file=paste(HOME,"/",lib,"/",lib,"_",type,"_rpm.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)

for (x in chr){
  file.remove(paste(HOME,"/",lib,"/",lib,".",x,"_",type,"_rpm.txt",sep=""))
}
