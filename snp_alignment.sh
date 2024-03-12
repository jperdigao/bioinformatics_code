
#arg1 list file
#arg2 path to directory containing .out files e.g. /path/to/directory/out/
#arg3 coverage e.g. 20
#arg4 fraction e.g. 0.75
#arg5 path to directory containing BAM and VCF files e.g. /path/to/directory/bam/

echo "List File: $1" >> parameters.txt
echo "Coverage Cutoff: $3" >> parameters.txt
echo "Proportion cutoff: $4" >> parameters.txt

sed "s~BAMPATH1~$5~g" ./step2_all_TB.pl > ./step2_all_TB1.pl
sed "s~BAMPATH1~$5~g" ./step3_all_TB.sh > ./step3_all_TB1.sh
sed "s~BAMPATH1~$5~g" ./step4_ALL_TB.sh > ./step4_ALL_TB1.sh

rm ./out_files/*.out
rm "$5"*.out

perl ./step2_all_TB1.pl $1
sh ./step3_all_TB1.sh
sh ./step4_ALL_TB1.sh
cp "$5"*.out ./out_files/

SNPtable_corrector.R SNP.global.0.2.list $1 $2 $3 $4
SNPtable_filter_Mtb_nores.R SNP.global.0.2.list_corrected "$1"_corrected


rm step2_all_TB1.pl step3_all_TB1.sh step4_ALL_TB1.sh
rm COV_FILES_LIST_ALL.out global_list_snp list_all_indels.txt list_all_snps.txt TOTAL*

