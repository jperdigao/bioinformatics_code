#!/usr/bin/perl

use strict;
use Getopt::Long;
use Carp;
use warnings;

#chr pos REf A C G T(covs) Rev-comp (1/2)
if (@ARGV !=4){
    print "#-------------------------------------------------------------\n";
    print "usage: calling_alleles.pl <In 1> <In 2> <In 3> <In 4>\n\n";
    print "\t<In 1> = File with  SNP file names..(absolute path!) \n";
    print "\t<In 2> = minimum allele cov.\n";
    print "\t<In 3> = % of total coverage eg. 0.2 for 20% \n";
    print "\t<In 4> = Multi call? [yes/no] [no=majority call; yes=Multi call]\n";
    #print "\t<In 6> = uniqueness file\n";
    print "#-------------------------------------------------------------\n";
    exit;
}
my ($file, $MIN_COV,$percentTotalcov, $call_method) = @ARGV;

#check inputs
if ($percentTotalcov > 1 or $percentTotalcov < 0){
    print "In-3 should be between 0 and 1\n"; exit;
}
if ($call_method ne "yes" && $call_method ne "no"){
    print $call_method , "\n";
      print "In-4 should be either yes or  no\n"; exit;
}

#my $MIN_COV = 10;

my %joinedFile;
my %ref_hash;
my %cov_hash;
my @order_keys;
my $key_counter =0;

#my %uniq_hash = %{hash_uniq(\$UNIQ_file)};
#Read file that contains file names... File names should be in absolute path..
open (FIL, $file);
while (<FIL>)
{
    chomp;
    my $fileName = $_;
    if($fileName eq ""){next;}
    open (FH, $fileName);
    
    my @arr = <FH>;
    if ($arr[0] =~/^Chromosome/) {
        shift @arr;
    }
    foreach (@arr)
    {
        my @sp = split /\s+/, $_;
        my ($chr, $pos, $ref, $cov_A, $cov_C, $cov_G, $cov_T, $rev) = @sp;
        my $tot =  $cov_A + $cov_C + $cov_G + $cov_T;
        my $key = $chr.":".$pos;
        if ($key_counter ==0){
            push @order_keys, $key;
        }
        if (! defined  $ref_hash{$key}){
            $ref_hash{$key} = $ref;
        }
        my $c="";
        #my $max ="";
        #my ($max, $max_pos)  = max(@covs);
        if ($call_method eq "yes"){
            
            my ($p) = $percentTotalcov*$tot;
            my $cutoff;
            if ($p>=$MIN_COV){
                $cutoff = $p;
            }else{
                $cutoff=$MIN_COV;
            }
            if( $tot == 0 ) {
		$c.="-"; 
	    } else {
            if ($cov_A > $cutoff){
                $c.="A";
            }if ($cov_C > $cutoff){
                $c.="C";
            }if ($cov_G > $cutoff){
                $c.="G";
            }if ($cov_T > $cutoff){
                $c.="T";
            }}
       	if ($c eq "") { $c.="N";} 
	}
        else{
            #$c = majorityCall (\@covs,\$MIN_COV);
		if( $tot == 0 ) {
                $c.="-";
            } else {           
	     if($cov_A >= $cov_C && $cov_A >= $cov_G &&  $cov_A >= $cov_T){
                $c.="A";
                #$max = $cov_A;
            }if($cov_C >= $cov_A && $cov_C >= $cov_G &&  $cov_C >= $cov_T){
                $c.="C";
                #$max = $cov_C;
            }if($cov_G >= $cov_A && $cov_G >= $cov_C &&  $cov_G >= $cov_T){
                $c.="G";
                #$max = $cov_G;
            }if($cov_T >= $cov_A && $cov_T >= $cov_C &&  $cov_T >= $cov_G){
                $c.="T";
                #$max = $cov_T;
            }
		}  
       	if ($c eq "") { $c.="N";} 
	}
        
        if (defined $joinedFile{$key}){
            $joinedFile{$key}.="%".$c;   
        }
        else{
            $joinedFile{$key}=$c;
        }
    }
   $key_counter+=1; 
}

foreach (@order_keys){
    chomp;
    my $key = $_;
    my @s1 = split /:/, $key;
    my @s2 = split /%/, $joinedFile{$_};
    #print $key, " ", $joinedFile{$_}, "\n";
    print $s1[0], " ", $s1[1] , " ", $ref_hash{$_}, " ";
    foreach (@s2){ print $_, " ";}
    print "\n";
}
exit;


