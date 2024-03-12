#!/usr/bin/perl

#verificar prefixo do vcf antes de executar na linha 47
use strict;
use warnings;

if (@ARGV!=1){
	print "USAGE: <bam_file_list> \n";
	exit;
}

my $path3="BAMPATH1";
my $ref="NC000962_3.fasta";
my $path5="BAMPATH1";
my $samtools="/run/media/myco1/Data/jperdigao/biosoft/samtools-0.1.18/samtools";
my $bcftools="/run/media/myco1/Data/jperdigao/biosoft/samtools-0.1.18/bcftools/bcftools";
my $vcfutils="/run/media/myco1/Data/jperdigao/biosoft/samtools-0.1.18/bcftools/vcfutils.pl";

#parameters for joining snps from coverage files...
my $UNIQ_cutoff=10;
my  $MIN_allele_cov=5;
my  $P_totalCall= 0.2;
my  $multiCall ="yes" ;
my $uniq_file = $path5."uniq.data";
my $thres = 30;

my ($bam_list) = @ARGV;
open BAM, $bam_list or die "can't open file $bam_list\n";

my %snp_hash;
my %indel_hash;
my %coverage_hash;

my @rm = `rm -f COV_FILES_LIST_ALL.out`;
my $cov_file_list = "COV_FILES_LIST_ALL.out";
my $cov_file_list1 = "COV_FILES_LIST_ALL.out";

open(out_file,">list_all_snps.txt");
open(out_file1,">list_all_indels.txt");

while (<BAM>)
{
	chomp;
	my $id = $_;
	my $prefix = $id;
	my $strain = $id;
	my $snp_file = $path3.$prefix.".concordant.eff.vcf";
	
	if (-e $snp_file)
	{
		print "opening ",$snp_file,"\n";
		open FIL, $snp_file or die "can't open file $snp_file\n";
		my $cov_file = $path3.$strain.".coverage";
		if (-e $cov_file){
		my @c = `echo $strain.coverage >> $cov_file_list`;	
		} 
		while (<FIL>)
		{
			next if 1 .. 26;
			chomp;	
			my @sp = split /\s+/, $_;
			 my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $names1) = ($sp[0], $sp[1], $sp[2], $sp[3], $sp[4], $sp[5], $sp[6], $sp[7], $sp[8], $sp[9]);
			my $key = $chr.":".$pos;  
			my $value = $strain.":".$alt; # print $value , "\n"; exit;
			 my $ref1 = substr $info, 0,2; 
			if ($ref1 ne "IN"  && $qual>=$thres) 
			{
				if (defined $snp_hash{$key}){
					$snp_hash{$key} = $snp_hash{$key}."%".$value;
				}
				else{
					$snp_hash{$key} = $value;
				}
			}	 		
			else{
					if (defined $indel_hash{$key}){
           	                    		$indel_hash{$key} = $indel_hash{$key}."%".$value; 
                        		}
                        		else{
                               			 $indel_hash{$key} = $value;
                      			}
			}
		}close(FIL);
	}
	else{
		print "couldn't find  file \n";
	}
}
#### --- Write out SNP hash ----

my @k_s = keys %snp_hash;
@k_s = sort(@k_s);
foreach (@k_s) { 
chomp; 
print out_file $_, " ", $snp_hash{$_}, "\n";
}

close(out_file);

#### --- Write out INDEL hash ----

my @l_s = keys %indel_hash;
@l_s = sort(@l_s);
foreach (@l_s) {
chomp;
print out_file1 $_, " ", $indel_hash{$_}, "\n";
}

close(out_file1);


#### --- Construct temp coverage files --- ####

open( FIL,$cov_file_list1);

while (<FIL>)
{
chomp; 
my $fileName = $_;
if($fileName eq "") {next;}
my $fileNameL = $path3.$fileName; 
my $fileNameO = $path3.$fileName.".out";
print "opening ",$fileNameL,"\n";
open (FH,$fileNameL);
open (OUT,">$fileNameO");
my @arr = <FH>;
if ($arr[0] =~/^Chromosome/) {
shift @arr;
}
print "writing ",$fileNameO,"\n";
foreach (@arr){
	my @s = split /\s+/, $_;
	#my ($chr, $pos, $ref, $cov_A, $cov_C, $cov_G, $cov_T, $rev) = @sp;
	#my $key1 = $chr.":".$pos;
	if( exists $snp_hash{$s[0].":".$s[1]} ) { print OUT $_; } 
}

close(FH);
close(OUT);


}


close(FIL);

exit; 

