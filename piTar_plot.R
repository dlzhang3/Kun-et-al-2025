library(Biostrings)
library(ggplot2)
library(reshape2)
library(scales)
library(ggpubr)
library(ggsci)
library(dplyr)
library(parallel)

comp_job <- "DZ.npp-14smRNA"
HOME <- "/gpfs/data/hlee-lab/"

path <- paste(HOME,"piTar/",comp_job,"/",sep="")
dir.create(path)
relaxed_path <- paste(path,"/relaxed_rules/",sep="")
stringent_path <- paste(path,"/stringent_rules/",sep="")
dir.create(relaxed_path)
dir.create(stringent_path)
flank <- 100

comp_key <- read.csv(paste(HOME,"piTar/comp_config/",comp_job,".csv",sep=""))

set_names <- comp_key$set_names
libs <- comp_key$libs
labs <- comp_key$labs
groups <- comp_key$groups
group_labs <- unique(comp_key$group_labs)

std <- function(x) sd(x)/sqrt(length(x))

normalize <- function(set_name,lib){
  annotation <- read.table(paste(HOME,"results/",set_name,"/",lib,"/libstat/annotations.txt",sep=""),sep="\t")
  nreads <- sum(annotation$v0)
  norm <- 10e6/nreads
  norm
}

norm_libs <- vector(length=length(libs))
for (x in 1:length(libs)){
  norm_libs[x] <- normalize(set_names[x],libs[x])
}

relaxed_WBGene <- readRDS(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/relaxed_WBGene",sep=""))
stringent_WBGene <- readRDS(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/stringent_WBGene",sep=""))

pos_vector <- as.character(seq(from=-(flank),to=flank))

read_rpm <- function(libs,set_name,sub_WBGene_list,relaxed_WBGene,stringent_WBGene){
  df_relaxed <- readRDS(file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_relaxed",sep=""))
  df_relaxed <- data.frame(do.call(rbind,df_relaxed))
  df_relaxed <- df_relaxed[which(relaxed_WBGene%in%sub_WBGene_list),]
  df_stringent <- readRDS(file=paste(HOME,"scratch/piRNA_sites_",set_name,"/",libs,"/",libs,"_stringent",sep=""))
  df_stringent <- data.frame(do.call(rbind,df_stringent))
  df_stringent <- df_stringent[which(stringent_WBGene%in%sub_WBGene_list),]
  
  
  ls <- list(df_relaxed,
             df_stringent)
  
  ls
  
}

build_result <- function(sub_WBGene_list,basename){
  
  relaxed_list <- list()
  stringent_list <- list()
  for (x in 1:length(libs)){
    sites <- read_rpm(libs=libs[x],
                      set_name=set_names[x],
                      sub_WBGene_list=sub_WBGene_list,
                      relaxed_WBGene,
                      stringent_WBGene)
    sites_relaxed <- sites[[1]]
    sites_stringent <- sites[[2]]
    
    relaxed_df <- data.frame(sites_relaxed)
    relaxed_sum <- data.frame(RPM=colSums(relaxed_df,na.rm=T)*norm_libs[x]/nrow(relaxed_df),
                              pos=pos_vector,
                              var=labs[x])
    relaxed_list[[x]] <- relaxed_sum
    
    stringent_df <- data.frame(sites_stringent)
    stringent_sum <- data.frame(RPM=colSums(stringent_df,na.rm=T)*norm_libs[x]/nrow(stringent_df),
                              pos=pos_vector,
                              var=labs[x])
    stringent_list[[x]] <- stringent_sum
    
  }
  
  refined_list <- list()
  for (i in unique(groups)){
    rpm_list <- list()
    index=1
    for (y in which(groups==i)){
      rpm_list[[index]] <- relaxed_list[[y]]$RPM
      index=index+1
    }
    rpm_df <- do.call(cbind,rpm_list)
    rpm_refined <- rowMeans(rpm_df,na.rm=T)
    sd_refined <- apply(rpm_df,1,sd,na.rm=T)
    
    refined_sum <- data.frame(RPM=rpm_refined,
                              SD=sd_refined,
                              pos=pos_vector,
                              var=group_labs[i])
    
    refined_list[[i]] <- refined_sum

  }
  
  relaxed_df <- do.call(rbind,refined_list)
  
  relaxed_df$SD[is.na(relaxed_df$SD)] <- 0
  
  relaxed_df$pos <- factor(relaxed_df$pos,
                           levels=pos_vector,
                           ordered=T)
  relaxed_df$var <- factor(relaxed_df$var,
                           levels=group_labs,
                           ordered=T)
  
  ymin <- min(relaxed_df$RPM)
  ymax <- ymin+(max(relaxed_df$RPM)-min(relaxed_df$RPM))*0.01
  
  gg_relaxed <- ggplot(data=relaxed_df,aes(x=pos,y=RPM,color=var,group=var))+
    theme_bw()+theme(aspect.ratio=1,
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     legend.position='top', 
                     legend.justification='left')+
    geom_line()+
    geom_ribbon(aes(x=pos,y=RPM,ymin=RPM-SD,ymax=RPM+SD,fill=var,group=var),na.rm=T,alpha=0.2,color=NA)+
    scale_color_brewer(palette="Set1")+
    scale_fill_brewer(palette="Set1")+
    scale_x_discrete(breaks=c("-50","-30","-10","10","30","50"))+
    labs(x="Position relative to 21U site",y="rpm/site")+
    coord_cartesian(xlim=c(((length(pos_vector)-1)/2-50),((length(pos_vector)-1)/2+52)))+
    annotate(geom="text",x=(length(pos_vector)-1)/2+12,y=ymax,label="U",color="darkgreen",size=3)
  
  rect_x <- 0
  rect_y <- rep(1,times=length(rect_x))
  rect_w <- 21
  
  df.rect <- data.frame(rect_x,
                        rect_y,
                        rect_w)
  
  piRNA <- ggplot(data=df.rect,aes(xmin=rect_x,xmax=rect_x+rect_w,ymin=rect_y,ymax=rect_y+1))+
    theme_void()+theme(legend.position = 0,
                       legend.title = element_blank())+
    scale_x_discrete(breaks=seq(from=df.rect$rect_x[1],to=df.rect$rect_x[nrow(df.rect)],by=100))+
    coord_cartesian(xlim=c(df.rect$rect_x[1],(df.rect$rect_x[nrow(df.rect)]+df.rect$rect_w[nrow(df.rect)])))+
    geom_rect(fill="darkgreen",color="darkgreen")
  piRNA.grob <- ggplotGrob(piRNA)
  
  
  
  gg_relaxed_full <- gg_relaxed+
    annotation_custom(
      grob=piRNA.grob,
      xmin=(length(pos_vector)-1)/2-9,
      xmax=(length(pos_vector)-1)/2+11,
      ymin=ymin,
      ymax=ymax
    )
  
  # ggsave(plot=gg_relaxed_full,
  #        file=paste("gg_relaxed_",basename,".png",sep=""),
  #        path=relaxed_path,
  #        dpi=300,
  #        height=5,
  #        width=5,
  #        device="png")
  
  ggsave(plot=gg_relaxed_full,
         file=paste("gg_relaxed_",basename,".pdf",sep=""),
         path=relaxed_path,
         dpi=300,
         height=5,
         width=5,
         device="pdf")
  
  
  refined_list <- list()
  for (i in unique(groups)){
    rpm_list <- list()
    index=1
    for (y in which(groups==i)){
      rpm_list[[index]] <- stringent_list[[y]]$RPM
      index=index+1
    }
    rpm_df <- do.call(cbind,rpm_list)
    rpm_refined <- rowMeans(rpm_df,na.rm=T)
    sd_refined <- apply(rpm_df,1,sd,na.rm=T)
    
    refined_sum <- data.frame(RPM=rpm_refined,
                              SD=sd_refined,
                              pos=pos_vector,
                              var=group_labs[i])
    
    refined_list[[i]] <- refined_sum
    
  }
  
  stringent_df <- do.call(rbind,refined_list)
  
  stringent_df$SD[is.na(stringent_df$SD)] <- 0
  
  stringent_df$pos <- factor(stringent_df$pos,
                             levels=pos_vector,
                             ordered=T)
  stringent_df$var <- factor(stringent_df$var,
                             levels=group_labs,
                             ordered=T)
  
  ymin <- min(stringent_df$RPM)
  ymax <- ymin+(max(stringent_df$RPM)-min(stringent_df$RPM))*0.01
  
  gg_stringent <- ggplot(data=stringent_df,aes(x=pos,y=RPM,color=var,group=var))+
    theme_bw()+theme(aspect.ratio=1,
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     legend.position='top', 
                     legend.justification='left')+
    geom_line()+
    geom_ribbon(aes(x=pos,y=RPM,ymin=RPM-SD,ymax=RPM+SD,fill=var,group=var),na.rm=T,alpha=0.2,color=NA)+
    scale_color_brewer(palette="Set1")+
    scale_fill_brewer(palette="Set1")+
    scale_x_discrete(breaks=c("-50","-30","-10","10","30","50"))+
    labs(x="Position relative to 21U site",y="rpm/site")+
    coord_cartesian(xlim=c(((length(pos_vector)-1)/2-50),((length(pos_vector)-1)/2+52)))+
    annotate(geom="text",x=(length(pos_vector)-1)/2+12,y=ymax,label="U",color="darkgreen",size=3)
  
  rect_x <- 0
  rect_y <- rep(1,times=length(rect_x))
  rect_w <- 21
  
  df.rect <- data.frame(rect_x,
                        rect_y,
                        rect_w)
  
  piRNA <- ggplot(data=df.rect,aes(xmin=rect_x,xmax=rect_x+rect_w,ymin=rect_y,ymax=rect_y+1))+
    theme_void()+theme(legend.position = 0,
                       legend.title = element_blank())+
    scale_x_discrete(breaks=seq(from=df.rect$rect_x[1],to=df.rect$rect_x[nrow(df.rect)],by=100))+
    coord_cartesian(xlim=c(df.rect$rect_x[1],(df.rect$rect_x[nrow(df.rect)]+df.rect$rect_w[nrow(df.rect)])))+
    geom_rect(fill="darkgreen",color="darkgreen")
  piRNA.grob <- ggplotGrob(piRNA)
  
  
  
  gg_stringent_full <- gg_stringent+
    annotation_custom(
      grob=piRNA.grob,
      xmin=(length(pos_vector)-1)/2-9,
      xmax=(length(pos_vector)-1)/2+11,
      ymin=ymin,
      ymax=ymax
    )
  
  # ggsave(plot=gg_stringent_full,
  #        file=paste("gg_stringent_",basename,".png",sep=""),
  #        path=stringent_path,
  #        dpi=300,
  #        height=5,
  #        width=5,
  #        device="png")
  
  ggsave(plot=gg_stringent_full,
         file=paste("gg_stringent_",basename,".pdf",sep=""),
         path=stringent_path,
         dpi=300,
         height=5,
         width=5,
         device="pdf")
  
  
}

gene.info <- read.table(paste0(HOME,"results/master.annotation.txt"),sep="\t",header=T)

all.WBGene <- gene.info$Gene.ID
csr.WBGene <- gene.info$Gene.ID[which(gene.info$CSR.target==T)]
wago.WBGene <- gene.info$Gene.ID[which(gene.info$WAGO1.target==T | gene.info$WAGO9.target==T)]
# activated.WBGene <- read.table("Y:/WJC_2021/20210708_figures/desilenced_wago_glh14.txt",sep="\n")
# activated.WBGene <- activated.WBGene$V1
# silenced.WBGene <- read.table("Y:/WJC_2021/20210708_figures/missilenced_csr_glh14.txt",sep="\n")
# silenced.WBGene <- silenced.WBGene$V1

# From 20210701_WJC_figures.Rmd, lines 4110-4111
# glh.activated.WBGene <- readRDS("Y:/piTar/adult_desilenced")
# meg.activated.WBGene <- readRDS("Y:/piTar/embryo_desilenced")

build_result(sub_WBGene_list=all.WBGene,
             basename="all")
build_result(sub_WBGene_list=csr.WBGene,
             basename="csr")
build_result(sub_WBGene_list=wago.WBGene,
             basename="wago")
# build_result(sub_WBGene_list=activated.WBGene,
#              basename="activated")
# build_result(sub_WBGene_list=silenced.WBGene,
#              basename="silenced")
# build_result(sub_WBGene_list=glh.activated.WBGene,
#              basename="glh_activated")
# build_result(sub_WBGene_list=meg.activated.WBGene,
#              basename="meg_activated")



# To build RDS object that I need for this script

# relaxed_bed <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.relaxed_rule.txt",sep=""),sep="\t",header=T)
# stringent_bed <- read.table(paste(HOME,"piTar/tmp.c_elegans.WS230.cds_transcripts/target_sites.stringent_rule.txt",sep=""),sep="\t",header=T)
# 
# extract_WBGene <- function(bed,bed_row){
#   unlist(strsplit(bed[bed_row,3],split="gene="))[2]
# }
# 
# relaxed_WBGene <- lapply(X=c(1:nrow(relaxed_bed)),FUN=extract_WBGene,bed=relaxed_bed)
# relaxed_WBGene <- do.call(rbind,relaxed_WBGene)
# relaxed_WBGene <- as.vector(relaxed_WBGene)
# 
# stringent_WBGene <- lapply(X=c(1:nrow(stringent_bed)),FUN=extract_WBGene,bed=stringent_bed)
# stringent_WBGene <- do.call(rbind,stringent_WBGene)
# stringent_WBGene <- as.vector(stringent_WBGene)
# 
# saveRDS(object=relaxed_WBGene,"Y:/piTar/tmp.c_elegans.WS230.cds_transcripts/relaxed_WBGene")
# saveRDS(object=stringent_WBGene,"Y:/piTar/tmp.c_elegans.WS230.cds_transcripts/stringent_WBGene")
