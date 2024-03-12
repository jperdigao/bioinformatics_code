## DECOMPRESS FASTQ FILES (PAIRED-END SEQUENCING) ##

gunzip FILE_1.fastq.gz

gunzip FILE_2.fastq.gz

## TRIM READS BY RUNNING TRIMMOMATIC IN PAIRED-END MODE ##

java -classpath trimmomatic.jar org.usadellab.trimmomatic.TrimmomaticPE -phred33 FILE_1.fastq FILE_2.fastq FILE_1_PAIRED_OUTPUT FILE_1_UNPAIRED_OUTPUT FILE_2_PAIRED_OUTPUT FILE_2_UNPAIRED_OUTPUT LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLIN:36

## INDEX REFERENCE GENOME ##

bwa index -a is NC000962_3.fasta

## MAP READS ##

bwa mem -M NC000962_3.fasta FILE_1_PAIRED_OUTPUT FILE_2_PAIRED_OUTPUT > FILE_PAIRED.sam

bwa mem -M NC000962_3.fasta FILE_1_UNPAIRED_OUTPUT > FILE_1_UNPAIRED.sam

bwa mem -M NC000962_3.fasta FILE_2_UNPAIRED_OUTPUT > FILE_2_UNPAIRED.sam

## CONVERT FILES FROM SAM FORMAT TO BAM FORMAT ##

java -jar picard.jar SamFormatConverter VALIDATION_STRINGENCY=LENIENT I=FILE_PAIRED.sam O=FILE_PAIRED.bam

java -jar picard.jar SamFormatConverter VALIDATION_STRINGENCY=LENIENT I=FILE_1_UNPAIRED.sam O=FILE_1_UNPAIRED.bam

java -jar picard.jar SamFormatConverter VALIDATION_STRINGENCY=LENIENT I=FILE_2_UNPAIRED.sam O=FILE_2_UNPAIRED.bam

## MERGE BAM FILES ##

java -jar picard.jar MergeSamFiles VALIDATION_STRINGENCY=LENIENT I=FILE_PAIRED.bam I=FILE_1_UNPAIRED.bam I=FILE_2_UNPAIRED.bam O=FILE.bam

## SORT BAM FILES ##

java -jar picard.jar SortSam VALIDATION_STRINGENCY=LENIENT I=FILE.bam O=FILE.bam SORT_ORDER=coordinate

## REMOVE DUPLICATE READS ##

java -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true I=FILE.bam O=FILE.bam METRICS_FILE=METRICS

java -jar picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=FILE.bam O=FILE.bam LB=FILE PL=illumina PU=FILE SM=FILE

## RE-SORT BAM FILES ##

java -jar picard.jar SortSam VALIDATION_STRINGENCY=LENIENT I=FILE.bam O=FILE.bam SORT_ORDER=coordinate

## INDEX BAM FILE ##

java -jar picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT I=FILE.bam

## REALIGN TARGETS AND INDELS ##

java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R NC000962_3.fasta -I FILE.bam -O TARGET_INTERVALS -allowPotentiallyMisencodedQuals #--fix_misencoded_quality_scores

java -jar GenomeAnalysisTK.jar -T IndelRealigner -R NC000962_3.fasta -I FILE.bam -targetIntervals TARGET_INTERVALS -O FILE.bam -allowPotentiallyMisencodedQuals #--fix_misencoded_quality_scores

## RE-INDEX BAM FILE ##

java -jar picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT I=FILE.bam

## VARIANT CALLING ##

java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R NC000962_3.fasta -I FILE.bam -ploidy 1 -glm BOTH -O FILE_GATK.vcf -allowPotentiallyMisencodedQuals

samtools mpileup -A -Q 23 -d 2000 -C 50 -ugf NC000962_3.fasta FILE.bam | bcftools call -O v -vm -o FILE_SAMTOOLS.vcf

## SELECT CONCORDANT VARIANTS ##

java -jar GenomeAnalysisTK.jar -T SelectVariants -R NC000962_3.fasta -V FILE_GATK.vcf --concordance FILE_SAMTOOLS.vcf -o FILE.vcf

## GENERATE A COVERAGE FILE CONTAINING COVERAGE ACROSS ALL POSITIONS IN THE REFERENCE GENOME

echo FILE > list.coverage

perl ./coverage.pl list.coverage 

rm list.coverage


############################################################################
## UPON CREATION OF BAM AND VCF FILES FOR ALL ISOLATES A SNP ALIGNEMENT CAN BE PRODUCED BY FILTER SHARED HIGH-QUALITY GENOMIC POSITIONS AND GENERATE AN ALIGNMENT ##
# The script takes a list of sample ids
# For a minimum depth coverage of 20 and allelic depth fraction of 0.75:

sh snp_alignment.sh list.txt /path/to/directory/out/ 20 0.75 /path/to/directory/bam/


############################################################################
## GENERATE PHYLOGENETIC TREE FROM THE FASTA FORMATED SNP ALIGNMENT##

iqtree -alrt 1000 -s ALIGNMENT.fasta
