#!/usr/bin/perl
unless(exists $ARGV[0])
{
	print "=======================================================\n";
	print "          Matrix to score (for paired predictions)\n";
	print "=======================================================\n";
	print " -mtx      raw gremlin matrix [required]\n";
	print " -div      length of protein A [required]\n";
	print " -cut      file mapping the full length sequence\n"; 
	print "           to the cut sequence [required]\n";
	print " -seq_len  number of sequences per length [required]\n";
	print "-------------------------------------------------------\n";
	print " -apcd     where to print APC corrected mtx [optional]\n"; 
	print "=======================================================\n";
	die();
}
sub mean
{
        my($data) = @_;
        if(not @$data){return 0;}
	else
	{
		my $total = 0;
		foreach (@$data){$total += $_;}
		my $average = $total / @$data;
		return $average;
	}
}
while ($arg = shift())
{
	if ($arg =~ m/^-/)
	{
		$arg = lc($arg);
		if ($arg eq "-mtx") {$mtx = shift(); next;}
		if ($arg eq "-div") {$div = shift(); next;}
		if ($arg eq "-cut") {$cut = shift(); next;}
		if ($arg eq "-seq_len") {$sl = shift(); next;}
		if ($arg eq "-apcd") {$apcd = shift(); next;}
	}
}
unless(-e $mtx){die("-mtx '$mtx' not found");}
unless(-e $cut){die("-cut '$cut' not found");}
unless(defined $div){die("-div '$div' not defined");}
unless(defined $sl){die("-seq_len '$sl' not defined");}
my $p_cut = 0.464276**($sl+0.995598) + 0.955303;
my $p_slope = 9.764891;

my @c2f;
my @c2i;
my @f2c;

my @FULL;my @CUT;
open(CUT,$cut);
my $n = 0;
while($line = <CUT>)
{
	chomp($line);
	if($n == 0){@FULL = split (//,$line);}
	if($n == 1){@CUT = split (//,$line);}
	$n++;
}
close(CUT);
my $f = 0;my $c = 0;
while(exists $FULL[$f])
{
	if($CUT[$f] ne "-")
	{
		$c2f[$c] = $f;
		$c2i[$c] = $CUT[$f];
		$c++;
	}
	$f++;
}

my @row_A_mean;
my @row_B_mean;
my @row_C_mean;
my @row_mean;
my @val_A;my @val_B;my @val_C;my @val;
my @MTX;
my @MODE;
my $i = 0;
open(MAT,$mtx);
while($line = <MAT>)
{
	chomp($line);
	my @tmp = split(/[\s\t]+/,$line);
	my @rows_A;
	my @rows_B;
	my @rows_C;
	my @rows;
	my $j = 0;
	foreach $t (@tmp)
	{
		$MTX[$i][$j] = $t;
		my $n = $i; if(exists $c2f[$i]){$n = $c2f[$i];}
		my $m = $j; if(exists $c2f[$j]){$m = $c2f[$j];}

		my $mode = "C";
		if($n < $div and $m < $div){$mode = "A";}
		elsif($n >= $div and $m >= $div){$mode = "B"}

		if($mode eq "A" and $i != $j){push(@rows_A,$MTX[$i][$j]);push(@val_A,$MTX[$i][$j]);}
		if($mode eq "B" and $i != $j){push(@rows_B,$MTX[$i][$j]);push(@val_B,$MTX[$i][$j]);}
		if($mode eq "C" and $i != $j){push(@rows_C,$MTX[$i][$j]);push(@val_C,$MTX[$i][$j]);}
		if($i != $j){push(@rows,$MTX[$i][$j]);push(@val,$MTX[$i][$j]);}

		$MODE[$i][$j] = $mode;
		$j++;
	}
	$row_A_mean[$i] = mean(\@rows_A);
	$row_B_mean[$i] = mean(\@rows_B);
	$row_C_mean[$i] = mean(\@rows_C);
	$row_mean[$i] = mean(\@rows);
	$i++;
}
close(MAT);

my $A_mean = mean(\@val_A);
my $B_mean = mean(\@val_B);
my $C_mean = mean(\@val_C);
my $mean = mean(\@val);

my @data;
my @val;
my $d = 0;
my $sc_max = 0;
my $i = 0;
while(exists $MTX[$i])
{
	my $j = $i + 1;
	while(exists $MTX[$i][$j])
	{
		my $apc = 1;
		if($MODE[$i][$j] eq "A"){$apc *= ($row_A_mean[$i]*$row_A_mean[$j])/$A_mean;}
		if($MODE[$i][$j] eq "B"){$apc *= ($row_B_mean[$i]*$row_B_mean[$j])/$B_mean;}
		if($MODE[$i][$j] eq "C"){$apc *= ($row_C_mean[$i]*$row_C_mean[$j])/$C_mean;}

		$MTX[$i][$j] = $MTX[$i][$j] - $apc;
		$MTX[$j][$i] = $MTX[$i][$j];

		if($MODE[$i][$j] eq "C")
		{
			if($MTX[$i][$j] > $sc_max){$sc_max = $MTX[$i][$j];}
			@{$data[$d]} = ($MTX[$i][$j],$MODE[$i][$j],$c2f[$i],$c2f[$j],$c2i[$i],$c2i[$j]);
			push(@val,$MTX[$i][$j]);
			$d++;
		}
		elsif($c2f[$j]-$c2f[$i] >= 3)
		{
			@{$data[$d]} = ($MTX[$i][$j],$MODE[$i][$j],$c2f[$i],$c2f[$j],$c2i[$i],$c2i[$j]);
			push(@val,$MTX[$i][$j]);
			$d++;
		}
		$j++;
	}
	$i++;
}

my $cut_len = int($i*1.5);

@val = sort {$b <=> $a} @val;
my @sub_val = @val[0..($cut_len-1)];

my $mean = mean(\@sub_val);
@data = sort { $b->[0] <=> $a->[0] } @data;

$sc_max = sqrt($sc_max/$mean);
my $p_max = 1/(1+exp(-(($sc_max-$p_cut)*$p_slope)));

print "i\tj\tgene\ti_id\tj_id\tr_sco\ts_sco\tp_sco\n";
my $n = 0;
while($n < $cut_len)
{
	my $mode = $data[$n][1];
	my $i = $data[$n][2]+1;
	my $j = $data[$n][3]+1;
	my $i_id = $i."_".$data[$n][4];
	my $j_id = $j."_".$data[$n][5];

	my $sco = sprintf("%.3f",$data[$n][0]);
	my $s_sco = sprintf("%.3f",$data[$n][0]/$mean);
	my $p_sco = "N/A";
	my $mode_ = $mode;
	if($mode eq "B")
	{
		$i_id = int($i-$div)."_".$data[$n][4];
		$j_id = int($j-$div)."_".$data[$n][5];
	}
	elsif($mode eq "C")
	{
		$p_sco = sprintf("%.3f",1/(1+exp(-(((sqrt($s_sco)*$p_max)-$p_cut)*$p_slope))));
		$i_id = $i."_".$data[$n][4];
		$j_id = int($j-$div)."_".$data[$n][5];
		$mode_ = "AB";
	}
	print $i."\t".$j."\t".$mode_."\t".$i_id."\t".$j_id."\t".$sco."\t".$s_sco."\t".$p_sco."\n";
	$n++;
}
if(defined $apcd)
{
	open(APCD,">$apcd");
	my $i = 0;
	while(exists $MTX[$i])
	{
		my $j = 0;
		while(exists $MTX[$i][$j])
		{
			my $val = 0;
			if($i != $j){$val = $MTX[$i][$j];}
			print APCD sprintf("%.6f",$val)." ";
			$j++;
		}
		print APCD "\n";
		$i++;
	}
	close(APCD);
}
