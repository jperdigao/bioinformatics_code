#!/bin/sh

path1=$PWD"/";
path2='BAMPATH1';

for isol in global;

do
echo $isol;
ls $path2*.coverage.out > "$isol"_list_snp;
perl calling_alleles.pl "$isol"_list_snp 10 0.2 no > gotemp
sed 's/[ \t]*$//' gotemp > SNP.$isol.0.2.list
done

rm gotemp
