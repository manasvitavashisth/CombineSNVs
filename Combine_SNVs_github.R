# Author: Manasvita Vashisth

rm(list=ls())
setwd("Set working directory")
library(tidyverse)
library(vcfR)
library(ggplot2)
library(data.table)
library(dplyr)

# Variable containing sample names
sample_names= c('Set sample names')

# Loop over each sample in list of sample names
for(i in 1:length(sample_names))
{
  
  # Initializes variables, setting each to NULL
  mutect2=NULL
  strelka=NULL
  varscan=NULL
  comb1=NULL
  comb2=NULL
  comb3=NULL
  # Reads the data from Varscan, Mutect2, and Strelka for each sample.
  varscan=as.data.table(fread(paste("Varscan path",sample_names[i],"/",sample_names[i],".snp.fpfiltered.hg38_multianno.txt",sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  mutect2=as.data.table(fread(paste("Mutect path",sample_names[i],"/pass_snvs.hg38_multianno.txt",sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  strelka=as.data.table(fread(paste("Strelka path",sample_names[i],"/results/variants/somatic.snvs.hg38_multianno.txt", sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE,na.strings=c(".", "NA")))
  
  # Renames columns of each dataset to make them consistent across all three variant callers.
  colnames(mutect2)[which(colnames(mutect2)=="Otherinfo1"):ncol(mutect2)]=c("Entry1","Entry2","Entry3","Chromosome_Number","Position","ID","Reference","Alternate","Quality","Mutect_Filter","Mutect_Info","Mutect_Format","Mutect_Normal","Mutect_Tumor")
  colnames(strelka)[which(colnames(strelka)=="Otherinfo1"):ncol(strelka)]=c("Entry1","Entry2","Entry3","Chromosome_Number","Position","ID","Reference","Alternate","Quality","Strelka_Filter","Strelka_Info","Strelka_Format","Strelka_Normal","Strelka_Tumor")
  colnames(varscan)[which(colnames(varscan)=="Otherinfo1"):ncol(varscan)]=c("Entry1","Entry2","Entry3","Chromosome_Number","Position","ID","Reference","Alternate","Quality","Varscan_Filter","Varscan_Info","Varscan_Format","Varscan_Normal","Varscan_Tumor")
  
  # Filters out variants that did not pass the specific variant caller's quality control (i.e. FILTER != PASS).
  mutect2=mutect2[mutect2$Mutect_Filter=="PASS",]
  strelka=strelka[strelka$Strelka_Filter=="PASS",]
  varscan=varscan[varscan$Varscan_Filter=="PASS",]
  
  # Creates a unique identifier for each variant, concatenating Chromosome, Start Position, Reference Allele, and Alternate Allele with '_'.
  mutect2$varID=paste0(mutect2$Chr,"_",mutect2$Start,"_",mutect2$Ref ,"_",mutect2$Alt)
  strelka$varID=paste0(strelka$Chr,"_",strelka$Start,"_",strelka$Ref,"_",strelka$Alt)
  varscan$varID=paste0(varscan$Chr,"_",varscan$Start,"_",varscan$Ref,"_",varscan$Alt)
  
  # Initialize columns in each dataframe to NA in place for the upcoming merge.
  mutect2[,c("Strelka_Filter","Strelka_Info","Strelka_Format","Strelka_Normal","Strelka_Tumor","Varscan_Filter","Varscan_Info","Varscan_Format","Varscan_Normal","Varscan_Tumor","S_GT","S_Depth","S_FDP","S_SDP","S_SUBDP","S_AU","S_CU","S_GU","S_TU"):=NA]
  strelka[,c("Mutect_Filter","Mutect_Info","Mutect_Format","Mutect_Normal","Mutect_Tumor","Varscan_Filter","Varscan_Info","Varscan_Format","Varscan_Normal","Varscan_Tumor"):=NA]
  varscan[,c("Strelka_Filter","Strelka_Info","Strelka_Format","Strelka_Normal","Strelka_Tumor","Mutect_Filter","Mutect_Info","Mutect_Format","Mutect_Normal","Mutect_Tumor","S_GT","S_Depth","S_FDP","S_SDP","S_SUBDP","S_AU","S_CU","S_GU","S_TU"):=NA]
  
  #Calculate Strelka VAF
  
  strelka <- separate(strelka, col = Strelka_Tumor, into = c("S_GT","S_Depth","S_FDP","S_SDP","S_SUBDP","S_AU","S_CU","S_GU","S_TU"), sep = ":", fill = "right")
  strelka_tier1=data.frame(matrix(nrow = nrow(strelka), ncol = 4))
  colnames(strelka_tier1)=c("A","C","G","T")
  strelka_tier1$A=as.numeric(str_extract(strelka$S_AU,"^([^,]+)"))
  strelka_tier1$C=as.numeric(str_extract(strelka$S_CU,"^([^,]+)"))
  strelka_tier1$G=as.numeric(str_extract(strelka$S_GU,"^([^,]+)"))
  strelka_tier1$T=as.numeric(str_extract(strelka$S_TU,"^([^,]+)"))
  a=c("A"=1,"C"=2,"G"=3,"T"=4)
  # Calculate total counts
  total_counts = rowSums(strelka_tier1)
  
  # Get the counts for the alternate allele
  alt_counts = strelka_tier1[cbind(1:nrow(strelka_tier1), match(strelka$Alt, names(a)))]
  
  # Calculate VAF
  strelka$S_VAF = alt_counts / total_counts
    comb1=bind_rows(mutect2[,c(1:129,143)],strelka[,c(1:129,151)],varscan[,c(1:129,143)])
  comb2=distinct(comb1,varID,.keep_all=TRUE)
  comb3=merge(comb2,mutect2[,138:143],by="varID",all=TRUE)
  comb3=merge(comb3,varscan[,138:143],by="varID",all=TRUE)
  comb3=merge(comb3,strelka[,c(138:151,162)],by="varID",all=TRUE)
  
  # Adds a flag/indicator for each variant caller if the variant has PASSed their respective FILTER
  comb3$mutect_flag=as.integer(comb3$varID %in% mutect2$varID)
  comb3$strelka_flag=as.integer(comb3$varID %in% strelka$varID)
  comb3$varscan_flag=as.integer(comb3$varID %in% varscan$varID)
  
  #Tissue extracting VAF and Read Depth
  df_split <- separate(comb3, col = Mutect_Tumor, into = c("MT_GT","MT_AD","MT_VAF","MT_Depth","MT_F1R2","MT_F2R1","MT_SB"), sep = ":", fill = "right")
  df_split <- separate(df_split, col = Varscan_Tumor, into = c("V_GT","V_GQ","V_Depth","V_RD","V_AD","V_VAF1","V_DP4"), sep = ":", fill = "right")
  df_split$V_VAF=as.numeric(str_extract(df_split$V_VAF1,"^([^%]+)"))/100
  df_split$MT_VAF=as.double(df_split$MT_VAF)
  df_split$MT_Depth=as.numeric(df_split$MT_Depth)
  df_split$S_Depth=as.numeric(df_split$S_Depth)
  df_split$V_Depth=as.numeric(df_split$V_Depth)
  df_split$V_RD=as.numeric(df_split$V_RD)
  df_split$V_AD=as.numeric(df_split$V_AD)
  df_split$S_FDP=as.numeric(df_split$S_FDP)
  df_split$S_SDP=as.numeric(df_split$S_SDP)
  df_split$S_SUBDP=as.numeric(df_split$S_SUBDP)
  df_split=df_split[(df_split$mutect_flag+df_split$strelka_flag+df_split$varscan_flag)>1,]
  df_split=df_split[df_split$gnomAD_genome_ALL<0.1 & df_split$ExAC_ALL<0.1 | is.na(df_split$gnomAD_genome_ALL) | is.na(df_split$ExAC_ALL),]
  df_split$sample_name=rep(sample_names[i],nrow(df_split))
  df_split$depth=rowMeans(df_split[,c('MT_Depth','S_Depth','V_Depth')],na.rm = TRUE)
  df_split$vaf=rowMeans(df_split[,c('V_VAF','S_VAF','MT_VAF')],na.rm=TRUE)
  df_split$patient=str_extract(df_split$sample_name, "[^_]+")
  df_split=df_split[df_split$depth>=20 & df_split$vaf>=0.1,]
  fwrite(df_split,file=paste(sample_names[i],"_Combined_TumorMutationCalls.txt",sep=""), sep="\t",col.names = TRUE)
}