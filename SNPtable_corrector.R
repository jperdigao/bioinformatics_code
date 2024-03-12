#!/usr/bin/env Rscript

cat("SNPtable_corrector v1.0 - - - Coverage-based correction of SNP Tables (manual pipeline)",sep="\n")
cat("",sep="\n")
cat("Author: João Perdigão",sep="\n")
cat("",sep="\n")
cat("",sep="\n")
cat("Usage: SNPtable_corrector.R [SNP tab file] [ordered strain list] [out files directory] [coverage cutoff] [fraction_cutoff]",sep="\n")
cat("",sep="\n")


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  cat("",sep="\n")
  cat("",sep="\n")
  cat("",sep="\n")
  cat("Usage: SNPtable_corrector.R [SNP tab file] [ordered strain list] [coverage cutoff] [fraction_cutoff]",sep="\n")
  cat("",sep="\n")
  stop("ERROR: This script requires a SNP table file...", call.=FALSE)
}



#####################

#setwd("D:\\lisbon_porto\\snp_table\\snp_table_caneti")

#cov_cutoff<-c(20)
#fraction_cutoff<-(0.75)

cov_cutoff<-as.numeric(args[4])
cat("",sep="\n")
cat(paste("Using coverage cutoff: ", args[4], " reads",sep=""))
cat("",sep="\n")

fraction_cutoff<-as.numeric(args[5])
cat("",sep="\n")
cat(paste("Using fraction cutoff: ", args[5],sep=""))
cat("",sep="\n")

cat("",sep="\n")
cat(paste("Using coverage out file directory: ", args[3],sep=""))
cat("",sep="\n")



snp_table_file<-args[1]
table1<-read.table(snp_table_file,sep=" ",header=F)

list<-read.table(args[2],sep=";",header=F)

cat("",sep="\n")
cat(paste("SNP table loaded:", args[1],sep=""))
cat("",sep="\n")

names(table1)[1:3]<-c("CHR","POS","REF")

names(table1)[4:length(names(table1))]<-as.vector(list$V1)

table2<-table1

strains_to_remove<-c()

for (strain in as.vector(list$V1)) {

  cat(paste("Looking at strain: ", strain,sep=""))

  coverage_file<-c(paste(args[3],strain,".coverage.out",sep=""))	
  
  strain_coverage<-read.table(coverage_file,sep="\t",header=F)
  
  names(strain_coverage)<-c("CHR","POS","REF","A","C","G","T","TOT")

  cat("",sep="\n")
  cat(paste("Strain file: ", coverage_file,sep=""))
  cat("",sep="\n")
  
  strain_coverage$TOT2<-strain_coverage$A+strain_coverage$C+strain_coverage$G+strain_coverage$T
  
  strain_coverage$fraction_tot<-strain_coverage$TOT2/strain_coverage$TOT
  
  #missed_positions<-as.vector(strain_coverage[strain_coverage$fraction_tot < 0.75,c("POS")])
  strain_coverage$max_cov<-apply(strain_coverage[,4:7], 1, function(x) max(x))
  strain_coverage$fraction_max<-strain_coverage$max_cov/strain_coverage$TOT
    
  missed_positions<-as.vector(na.omit(strain_coverage[strain_coverage$fraction_max < as.numeric(fraction_cutoff),c("POS")]))
  
  table2[,c(strain)]<-as.character(table2[,c(strain)])
  
  table2[table2$POS %in% strain_coverage$POS[strain_coverage$TOT < as.numeric(cov_cutoff)],c(strain)]<-c("-")

  table2[table2$POS %in% missed_positions,c(strain)]<-c("-")	

  #for (i in as.vector(table2$POS)) {
    
  #  if (strain_coverage[strain_coverage$POS == c(i),c("TOT")] < cov_cutoff) {table2[table2$POS == c(i),c(strain)]<-c("-")}
  
  #  if (i %in% missed_positions) {table2[table2$POS == c(i),c(strain)]<-c("-")}
    
    
    
    
    
  #}
  
  total_SNPs<- length(as.vector(table2$POS))
  missed_length<-length(as.vector(table2[table2[[strain]] == c("-"),c("POS")]))
  missed_fraction<-missed_length/total_SNPs
  
  if (missed_fraction > 0.1) {strains_to_remove<-c(strains_to_remove,strain)}
  
  
  cat("...done!")  
  cat("",sep="\n")  
}
    
no_strains<-length(names(table2))-3

#table3<-table2

#for (i in as.vector(table3$POS)) {
#  
#  pos_alleles<-as.vector(t(table3[table3$POS == c(i),-c(1,2,3)]))
#  
#  missed_calls_pos<-pos_alleles[pos_alleles == "-"]
#  
#  missed_pos_fraction<-length(missed_calls_pos)/no_strains
#  
#  if (missed_pos_fraction > 0.1) {table3<-table3[-(table3$POS == c(i)),]}
#
#  cat(paste("Looking missed calls at position: ", i,sep=""))
#  cat(paste("...", match(i,table3$POS),"/",length(as.vector(table3$POS))))  
#  cat("",sep="\n")
#}

################# Remove positions in excess of 10% missed calls 
table_ref<-table2

table_ref$miss_count<-rowSums(table_ref[,4:(length(names(table_ref))-1)] == "-")

table_ref$miss_count_ratio<-table_ref$miss_count/length(names(table_ref)[3:(length(names(table_ref))-1)])


table3<-table_ref[table_ref$miss_count_ratio <= 0.1,c(1:(length(names(table_ref))-2))]

cat(paste(length(table_ref$miss_count_ratio[table_ref$miss_count_ratio > 0.1])," Positions removed with an excess missed calls of  10%",sep=""))
cat("",sep="\n")

################

table4<-table3

table4[strains_to_remove]<-NULL

list_final<-data.frame(V1=names(table4)[4:length(names(table4))])

write.table(table4,paste(snp_table_file,"_corrected",sep=""),sep=" ",col.names = T,row.names = F,quote = F)

write.table(table3,paste(snp_table_file,"_corrected_allstrains",sep=""),sep=" ",col.names = T,row.names = F,quote = F)
  
write.table(list_final,paste(args[2],"_corrected",sep=""),sep=" ",col.names = F,row.names = F,quote = F)

cat(paste("Strains removed: ", strains_to_remove,sep="\n"))
