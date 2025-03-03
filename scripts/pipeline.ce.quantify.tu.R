#!$HOME/bin/Rscript

HOME <- getwd()

args<-commandArgs(TRUE)
intersect<-read.table(args[1],sep="\t")
out     <-args[2]
norm    <-as.numeric(args[3])
id.class<-args[4]

gene.id <- read.table(paste0(HOME,"/geneIDs.WS230"), sep="," );
colnames(gene.id)<-c("WBGene","name","cosmid")

intersect.sense<-intersect[intersect[,6]==intersect[,12],]
gene.sense.ppm<-as.data.frame(tapply(intersect.sense[,5],intersect.sense[,10],sum)*1000000/norm)
colnames(gene.sense.ppm)<-"ppm"
gene.sense.ppm<-merge(gene.sense.ppm,gene.id,by.x=0,by.y=id.class,all.x=T)
colnames(gene.sense.ppm)[0]<-id.class
write.table(gene.sense.ppm[,c(3,4,1,2)],paste(out,"sense","ppm",sep="."),sep="\t",quote=F,row.names=F,na="0")
write.table(intersect.sense[,c(4,5)],paste(out,"sense","ntm",sep="."),sep="\t",quote=F,row.names=F,col.names=F,na="0")

intersect.anti<-intersect[intersect[,6]!=intersect[,12],]
gene.anti.ppm<-as.data.frame(tapply(intersect.anti[,5],intersect.anti[,10],sum)*1000000/norm)
colnames(gene.anti.ppm)<-"ppm"
gene.anti.ppm<-merge(gene.anti.ppm,gene.id,by.x=0,by.y=id.class,all.x=T)
colnames(gene.anti.ppm)[0]<-id.class
write.table(gene.anti.ppm[,c(3,4,1,2)],paste(out,"anti","ppm",sep="."),sep="\t",quote=F,row.names=F,na="0")
write.table(intersect.anti[,c(4,5)],paste(out,"anti","ntm",sep="."),sep="\t",quote=F,row.names=F,col.names=F,na="0")
