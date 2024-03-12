#!/bin/sh

path1='BAMPATH1';
path2="$PWD""/"
filen='COV_FILES_LIST_ALL.out'
files=`awk '{ print $1 }' $path2$filen`;
echo $files;
files1=`head -1 $path2$filen`;
echo $files1;

awk '{print $1,$2,$3}' $path1$files1.out > TOTAL_COV_FILES_LIST_ALL; 

for file in $files

do

# awk '{print $1,$2,$3,$4+$5+$6+$7}' $path1$file.out > $path1$file.total;
awk '{print $4+$5+$6+$7}' $path1$file.out > temp2;
paste TOTAL_COV_FILES_LIST_ALL temp2 -d " " >> coverage1.out;  
mv coverage1.out TOTAL_COV_FILES_LIST_ALL 
rm temp2;
done

awk '{sum=0; for(i=4; i<=NF; i++){sum+=$i}; print $1,$2,$3,sum}'  TOTAL_COV_FILES_LIST_ALL > TOTAL_COV_ALL.txt
