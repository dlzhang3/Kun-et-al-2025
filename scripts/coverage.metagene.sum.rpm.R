#!$HOME/bin/Rscript

library(parallel)
library(data.table)

args <- commandArgs(T)

lib <- args[1]
type <- args[2]
chr <- args[3]

HOME <- getwd()

df.bedgraph <- fread(paste(HOME,"/",lib,"/",lib,".inserts.v0.",type,".norm.anti.depth.",chr,sep=""),sep="\t",header=F)

transcriptome <- read.table(paste(HOME,"/ce_WS230.coding_transcript.exon.merge.gene.bed",sep=""),sep="\t",header=F)

gene.list <- as.character(transcriptome[which(transcriptome[,1]==chr),4])

bin_cds_rpm <- function(depth.graph,WBGene,n.bin,transcriptome){
  row <- transcriptome[which(transcriptome[,4]==WBGene),]
  range <- vector(length=0L)
  exon.starts <- as.numeric(as.vector(unlist(strsplit(x=as.character(row[1,12]),split=","))))
  exon.lengths <- as.numeric(as.vector(unlist(strsplit(x=as.character(row[1,11]),split=","))))
  n.exons <- row[1,10]
  strand <- as.character(row[1,6])
  if (strand=="+"){
    cds.start <- row[1,2]
    for(ii in 1:n.exons){
      range <- c(range,(cds.start+exon.starts[ii]):(cds.start+exon.starts[ii]+exon.lengths[ii]))
    }
  } else if (strand=="-"){
    cds.start <- row[1,3]
    for(ii in 1:n.exons){
      range <- c(range,(cds.start-exon.starts[ii]):(cds.start-exon.starts[ii]-exon.lengths[ii]))
    }
  }

  
  if (length(range)>n.bin){
    bin.size <- length(range)/n.bin
    rel.vector <- list()
    
    rel.vector[[1]] <- range[1:bin.size]
    
    for (x in 2:n.bin){
      last <- which(range==rel.vector[[x-1]][bin.size])+1
      rel.vector[[x]] <- range[last:(last+bin.size-1)]
    }
    
    bin.values <- vector(length=n.bin)
    for (i in 1:n.bin){
      bin.values[i] <- sum(depth.graph[rel.vector[[i]],3])
    }
    total.length <- sum(bin.values)
    
    if (total.length>=5){
      df <- data.frame(bin.values)
      colnames(df) <- WBGene
      
      df
    }
  }
}

ls <- mclapply(FUN=bin_cds_rpm,X=gene.list,depth.graph=df.bedgraph,n.bin=100,transcriptome=transcriptome)

df <- data.frame(NULL)
for (x in 1:length(ls)){
  if (length(ls[[x]])>0){
    df <- rbind(df,t(ls[[x]]))
  }
}

write.table(x=df,file=paste(HOME,"/",lib,"/",lib,".",chr,"_",type,"_rpm.txt",sep=""),sep="\t",quote=F,col.names=F)
