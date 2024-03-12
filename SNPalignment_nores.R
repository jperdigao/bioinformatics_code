#!/usr/bin/env Rscript

cat("SNPtable_filter_Mtb v2.0 - - - From SNP Table manual pipeline",sep="\n")
cat("",sep="\n")
cat("Author: João Perdigão",sep="\n")
cat("",sep="\n")
cat("Usage: SNPtable_filter_Mtb.R [SNP tab file] [ordered strain list]",sep="\n")
cat("",sep="\n")

require(ape)
require(seqinr)
#require(phangorn)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("ERROR: This script requires a SNP table file...", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"  ##Include standard output file
}

cat("Reading M. tuberculosis library files...",sep="\n")
cat("",sep="\n")

ppe_positions<-read.csv("ppe_positions_all.txt",sep="\t",header=F)
database_uniq<-read.table("kmer.50.Reference.pre.txt",header=F,sep=" ")
res_pos<-read.csv("res_pos_extended_list.txt",sep="\t",header=T)


cat("Reading M. tuberculosis library files...OK!",sep="\n")
cat("",sep="\n")

names(database_uniq)<-c("POS","V2")
pos_rem<-database_uniq[database_uniq$V2 < 49,]


snp_table_file<-args[1]
snp_table<-read.csv(args[1],sep=" ",header=T)

names(snp_table)[c(1,2,3)]<-c("CHR","POS","REF")
len<-length(names(snp_table))
#strains<-read.csv(args[2],sep=";",header=F)
#names(snp_table)[c(4:len)]<-as.vector(strains$V1)
#cat("checkpoint 1 ok")
snp_table<-subset(snp_table, select=-c(CHR))


snp_table1<-snp_table[!snp_table$POS %in% ppe_positions$V1,]
cat("PE/PPE Positions removed...",sep="\n")
cat("",sep="\n")
snp_table2<-snp_table1[!snp_table1$POS %in% pos_rem$POS,]
cat("Positions with mapping quality below k-mer 49/50 removed...",sep="\n")
cat("",sep="\n")
snp_table2<-snp_table2[!snp_table2$POS %in% res_pos$pos,]
cat("Positions within resistance genes removed...",sep="\n")



cat("Converting to multifasta format...",sep="\n")
cat("",sep="\n")

snp_table3<-snp_table2 #[,-1]
#to remove "-" - coreSNPs (shared between all isolates)
#snp_table3<-as.data.frame(lapply(snp_table3,function(x) gsub("-",NA,x)))
snp_table3<-snp_table3[complete.cases(snp_table3),]

#to remove double alleles
for (i in 2:length((names(snp_table3)))) {
  
  for (a in 1:length(snp_table3$REF)){
    
    if (nchar(as.character(snp_table3[a,i])) > 1) {snp_table3[a,i]<-NA}
      }
}

snp_table3<-snp_table3[complete.cases(snp_table3),]

core_sites_length<-length(snp_table3$REF)

cat(paste("Total no. of core SNP sites...",core_sites_length,sep=""))
cat("",sep="\n")

#to write a csv SNP table filtered:
write.table(snp_table3,"SNP_table_filtered_nores.csv",sep=";",col.names=T,row.names=F,quote=F)
###
snp_table3<-snp_table3[,-1]

sequence<-t(sapply(strsplit(as.character(snp_table3[,1]),""),tolower))
#sequence<-t(strsplit(as.character(snp_table3[,1]),""))

row.names(sequence)<-names(snp_table3[1])
alignment<-as.DNAbin(sequence)

for (i in 2:length((names(snp_table3)))) {
  
  
  sequence2<-t(sapply(strsplit(as.character(snp_table3[,i]),""),tolower))        
  #row.names(sequence2)<-names(snp_table[i])
  
  alignment2<-as.DNAbin(sequence2)
  
  alignment<-rbind(alignment,alignment2)
  
  
}

row.names(alignment)<-names(snp_table3)
#row.names(alignment)[1]<-c("H37Rv") # adjust reference name

cat("Conversion to multifasta format... OK!",sep="\n")
cat("",sep="\n")
write.dna(alignment, "wgSNP_alignment_filtered_nores.fas", format = "fasta",colsep = "")

cat("Output FASTA file: wgSNP_alignment_filtered.fas",sep="\n")
