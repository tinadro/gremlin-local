#!/usr/bin/perl
unless(exists $ARGV[0])
{
	print "========================================================================\n";
	print "           remove gappy positions from alignment\n";
	print "========================================================================\n";
	print " -i        a3m/fas [required]\n";
	print " -percent  minimum % of non_gap characters per position [default=25]\n";
	print " -rm       positions to exclude [example: 0-30,40-100]\n";	
	print " -prefix   prefix added to files made [default=i]\n";	
	print "========================================================================\n";
	die();
}
my @cut;
while ($arg = shift())
{
	if ($arg =~ m/^-/)
	{
		if ($arg eq "-i") {$i = shift(); next;}
		if ($arg eq "-percent") {$percent = shift(); next;}
		if ($arg eq "-rm"){@cut = split(",",shift());next;}
		if ($arg eq "-prefix"){$prefix = shift();next;}
	}
}
my @CUT; my $n = 0;for my $c (@cut){@{$CUT[$n]} = split("-",$c);$n++;}

unless(-e $i){die("-i '$i' not found");}
unless(defined $prefix){$prefix = substr($i,0,-4);}
unless(defined $percent){$percent = 0.25;}
if($percent > 1){$percent = $percent/100;}

# read fasta
my @col;
my $n = 0;
my @seq;
my $s = 0;
open(A3M,$i);
while($line = <A3M>)
{
	chomp($line);
	if(length($line) > 0)
	{
		if(substr($line,0,1) eq ">")
		{
			if(exists $seq[$s]){$s++;}
			$seq[$s][0] = $line;
		}
		else
		{
			$line =~ s/[a-z]//g; #make fasta
			$seq[$s][1] .= $line;
		}
	}
}
close(A3M);
my $seq_ = $s;

my $len_ = length($seq[0][1]);
my $n = 0;
while(exists $seq[$n])
{
	my $m = 0; for my $p (split(//,$seq[$n][1])){$col[$m] .= $p; $m++;}
	$n++;
}
print "i\tnon_gap\tseq_len\n";
my @col_check;
my $m = 0;
while(exists $col[$m])
{

	my $gap = int($col[$m] =~ tr/-//);
	my $non_gap = ($seq_-$gap);

	if($non_gap/$seq_ >= $percent){$col_check[$m] = 1;}
	else{$col_check[$m] = 0;}

	my $n = 0;
	while(exists $CUT[$n])
	{
		if($m >= $CUT[$n][0] and $m <= $CUT[$n][1]){$col_check[$m] = 0;}
		$n++;
	}
	if($col_check[$m] == 0){$cut_seq .= "-";}
	else{$cut_seq .= substr($seq[0][1],$m,1);}
	print $m."\t".$col_check[$m]."\t".sprintf("%.3f",$non_gap/$seq_)."\t".sprintf("%.3f",$non_gap/$len_)."\n";
	$m++;
}
open(CUT,">$prefix.cut");
print CUT $seq[0][1]."\n".$cut_seq."\n";
close(CUT);
$cut_seq =~ s/-//g;
my $cut_len = length($cut_seq);
open(FASTA,">$prefix.cut.fas");
open(MSA,">$prefix.cut.msa");
my $n = 0;
my $N = 0;
while(exists $seq[$n])
{
	my $m = 0;
	my $new_seq;
	while($m < $len_)
	{
		if($col_check[$m] == 1)
		{
			my $p = substr($seq[$n][1],$m,1);
			$new_seq .= $p;
		}
		$m++;
	}
	print FASTA $seq[$n][0]."\n".$new_seq."\n";
	print MSA "S_".sprintf("%06d", ($N+1))."     ".$new_seq."\n";
	$N++;
	$n++;
}
close(FASTA);
close(MSA);
print STDERR "seq_len (number of sequences per length): ".sprintf("%.3f",$N/$cut_len)."\n";
