library(Biostrings)
library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)
library(ggsci)
library(dplyr)
library(parallel)

args <- commandArgs(T)

libs <- args[1]
set_name <- args[2]
rules_path <- as.numeric(as.character(args[3]))

HOME <- "/gpfs/data/hlee-lab/"

if (rules_path<51){
  rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.relaxed_rule_",rules_path,".txt",sep=""),sep="\t",header=F)
  rule_name <- paste("relaxed_",rules_path,sep="")
} else if (rules_path==51){
  rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.stringent_rule.txt",sep=""),sep="\t",header=T)
  rule_name <- "stringent"
} else if (rules_path==52){
  rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.common_stringent.txt",sep=""),sep="\t",header=T)
  rule_name <- "common_stringent"
} else if (rules_path==53){
  rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.common_relaxed.txt",sep=""),sep="\t",header=T)
  rule_name <- "common_relaxed"
} else if (rules_path>53){
  index <- rules_path-53
  rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/CLASH/target_sites.CLASH_",index,".txt",sep=""),sep="\t",header=T)
  rule_name <- paste("CLASH_",index,sep="")
}

# libs <- "WJC.20180806_RRS"
# set_name <- "WJC.20180806"

bedgraph_to_depth <- function(df.sub,range){
  depth <- rep(x=0,times=length(range))
  
  if (max(which(between(df.sub$V2,min(range),max(range))),0)>0){
    select <- df.sub[which(between(df.sub$V2,min(range),max(range))),]
    
    
    
    for (i in 1:nrow(select)){
      
      if (max(which(range==select[i,3]),0)==0){
        select[i,3] <- max(range)
      }
      
      depth[which(range==select[i,2]):which(range==select[i,3])] <- select[i,4]
    }
  }
  
  
  
  depth
}

# relaxed_rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.relaxed_rule.txt",sep=""),sep="\t",header=T)
# stringent_rules <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.stringent_rule.txt",sep=""),sep="\t",header=T)

get_rpm <- function(rules,rule_rownum,track,flank){
  site <- rules[rule_rownum,]
  gene <- as.character(site[,3])
  transcript_start <- (site[,4])-(flank-11)
  transcript_end <- site[,5]+(flank-10)
  
  range <- transcript_start:transcript_end
  
  df.sub <- track[which(track[,1]==gene),]
  
  depth <- bedgraph_to_depth(df.sub,range)
  
  # print(rule_rownum)
  
  depth
  
}


track <- read.table(paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,".inserts.v0.22G.norm.anti.bedgraph",sep=""),sep="\t",header=F)


ls <- mclapply(X=c(1:nrow(rules)),FUN=get_rpm,
               rules=rules,
               flank=100,
               track=track)

saveRDS(object=ls,
        file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_",rule_name,sep=""))

df <- do.call(rbind,ls)

write.table(x=df,
            file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_",rule_name,".txt",sep=""),
            quote=F,
            row.names=F,
            col.names=F)


# relaxed_ls <- mclapply(X=c(1:nrow(relaxed_rules)),FUN=get_rpm,
#                            rules=relaxed_rules,
#                            flank=50,
#                            track=track)
# 
# saveRDS(object=relaxed_ls,
#         file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_relaxed",sep=""))
# 
# relaxed_df <- do.call(rbind,relaxed_ls)
# 
# write.table(x=relaxed_df,
#             file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_relaxed.txt",sep=""),
#             quote=F,
#             row.names=F,
#             col.names=F)
# 
# 
# stringent_ls <- mclapply(X=c(1:nrow(stringent_rules)),FUN=get_rpm,
#                      rules=stringent_rules,
#                      flank=50,
#                      track=track)
# 
# saveRDS(object=stringent_ls,
#         file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_stringent",sep=""))
# 
# stringent_df <- do.call(rbind,stringent_ls)
# 
# write.table(x=stringent_df,
#             file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_stringent.txt",sep=""),
#             quote=F,
#             row.names=F,
#             col.names=F)
