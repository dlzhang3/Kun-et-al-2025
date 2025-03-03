#!$HOME/bin/Rscript

HOME <- getwd()

args<-commandArgs(TRUE)
type <- args[length(args)]

master <- read.table(paste(HOME,"/reference/master.",type,".txt",sep=""),sep="\t",header=T)

lib.data <- list()
for (i in 1:(length(args)-1)){
  lib.data[[i]]=read.table(paste(HOME,"/results/",args[i],"/rpm/",args[i],".inserts.xk.uniq.reads.",
                                 type,".ppm",sep=""),sep="\t",header=T)
}


for (x in 1:(length(args)-1)){
  master <- merge(master,lib.data[[x]],by.x="Gene.ID",by.y="Row.names",all.x=T,all.y=T)
  master <- master[,c(1:(3+x-1),((3+x-1)+3))]
}

colnames(master) <- c("Gene.ID","name","cosmid",
                                            args[1:(length(args)-2)])

write.table(master,paste("master.",type,".txt",sep=""),sep="\t",quote=F,row.names=F,na="0")
