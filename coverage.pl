#!/usr/bin/perl

use strict;
use warnings;
##use Carp;
##use POSIX qw(ceil floor);
## use Settings ;
use Getopt::Long ;

my $path= "./";
my $path2="./";
my $path1="./";
my $reference=$path1."NC000962_3.fasta";
my $samtools="samtools";
my $bcftools="bcftools";
my $vcfutils="vcfutils.pl";

if (@ARGV !=1){
    print "<run_lane file>\n";
    exit;
}

my ($file, $q, $mem) =@ARGV;
open (FH, $file) or "die can't open file\n";

my %hash;
my @s;
while (<FH>){
    chomp;
    my @sp = split /\s+/, $_;
    my ($sample) = @sp;
    my $dir = $sample;
    my $b = $sample.".bam";
   if (defined $hash{$sample}){
    $hash{$sample}.="%".$b;
   }
   else{
    $hash{$sample} = $b;
    push @s, $sample;
   }
}close (FH);

foreach (@s){
    chomp;
    my @sp = split /%/, $hash{$_};
    my $bam_files ="";
    my $n = scalar(@sp);
    foreach (@sp){
        chomp;
       if ($bam_files eq ""){ $bam_files= "$path$_";}
       else{
        $bam_files.=" $path$_";
       }
    }
    my $sample = $_;
    # print "$_  $n $bam_files\n";
    my $output = $path2.$sample.".coverage";
    print "$sample \n";
    my $cs1 ="$samtools mpileup -f $reference $bam_files > output.$sample";
    my $cs2 = "perl rawPileup_2_coverage.pl  output.$sample  $reference > $output";
    my @do = `$cs1 && $cs2 && rm -f output.$sample`;
    # my @do = `$cs1 && $cs2`;
}
exit;
