#!/software/bin/perl
use strict;
use warnings;
use Carp;
use POSIX qw(ceil floor);
use Getopt::Long;
my ($indels);
GetOptions ('indels' => \$indels);

if (@ARGV <2){
    print "<raw pileup file> <ref_file> [--indels]\n";
    exit;
}

if ($indels) {
    print STDERR "also printing indels\n";
}
my $p=$ARGV[0];
my $ref= $ARGV[1];
my ($r1, $r2) = hash_contigs($ref);
my %ref_hash = %{$r1};
my @chr_ids = @{$r2};

print "Chromosome\tPosition\tReference\tA\tC\tG\tT\tTotal";
print "\tInsertions\tDeletions" if $indels ;
print "\n" ;

foreach (@chr_ids)
{
    chomp;
    my %cov_hash;
    my $chr = $_;
    my $com = "awk \'{if (\$1==\"$chr\"){print }}\' $p";
    my @FH = `$com`;
    foreach (@FH)
    {
	chomp;
	my @sp = split /\s+/, $_;
	my $n = scalar (@sp);
	if ($n < 4 ){next;}
	my $chr = $sp[0];
	my $pos = $sp[1];
	my $ref = $sp[2];
	$ref = uc($ref);
	my $bases = $sp[$n-2];
	$bases =~s/[.,]/$ref/g;

	my $count = $sp[$n-3];
	my $insertions =0;
	my $deletions =0;
	if ($indels){
	    while ($bases =~/\+(\d+)[ACGT]/g){
		$insertions+=$1;
	    }
	    while ($bases=~/\-(\d+)[ACGT]/g){
		$deletions+=1;
	    }
	}
	$bases=~ s/[\+\-]\d+[ACGTN]//g;
	#print "$chr $pos $count $pil\n";
	my ($a, $c, $g, $t)=(0,0,0,0,0);
	$a = $bases =~ tr/A//;
	$c = $bases =~ tr/C//;
	$g = $bases =~ tr/G//;
	$t = $bases =~ tr/T//;
	my $sum = $a + $c+$g+$t;
	my $id = $chr."%".$pos;
	if ($indels){
	    $cov_hash{$id}="$chr\t$pos\t$ref\t$a\t$c\t$g\t$t\t$sum\t$insertions\t$deletions\n";
	}
	else{
	    $cov_hash{$id}="$chr\t$pos\t$ref\t$a\t$c\t$g\t$t\t$sum\n";
	}
    }
    my $len = length($ref_hash{$_});
    my @seq= split //,$ref_hash{$chr};
    for (my $i=1; $i<=$len; $i+=1){
	my $id = $chr."%".$i;
	if (defined $cov_hash{$id}){
	    print $cov_hash{$id};
	}
	else{
	    if (defined $seq[$i-1]){
		my $ref=uc($seq[$i-1]);
		if ($indels)
		{
		    print "$chr\t$i\t$ref\t0\t0\t0\t0\t0\t0\t0\n";
		}
		else{
		    print "$chr\t$i\t$ref\t0\t0\t0\t0\t0\n";
		}
	    }
	}
    }
   
}
exit;

#------------------------------------------------------
#---------------  hash_contigs
#------------------------------------------------------
sub hash_contigs{
  if(  @_ != 1 )
	{
	  print "Usage: hash_contigs contigs_file";
	  exit;
	}
  my $contigs_file = shift;
  #print "#	hashing contigs from the file.. $contigs_file\n";
  if( $contigs_file =~ /\.gz$/ ){
	  open(CONTIGS, $contigs_file, "gunzip -c $contigs_file |" ) or die "Cant open contigs file: $!\n";
	}else{
	  open( CONTIGS, $contigs_file) or die "Cant open contigs file: $!\n";
	}
  my %contigs_hash; # hash to store contig names and sequences
  my $contigName;
  my @ids;
  while (<CONTIGS>){
	  if (/^>(\S+)/) {
		$contigName=$1;
                push @ids, $contigName;
	} else {
		chomp;
		$contigs_hash{$contigName} .= $_;
	  }
	}
  close(CONTIGS);
  return (\%contigs_hash, \@ids);
}

sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

