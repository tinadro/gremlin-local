#!/usr/bin/perl
my $s = 0;
my $ma_nam;
my @order;
my %STO;
open(STO,$ARGV[0]);
while($line = <STO>)
{
	chomp($line);
	if(substr($line,0,1) ne "#")
	{
		my ($nam,$seq) = split(/[\t\s]+/,$line);
		if(defined $seq)
		{
			unless(exists $STO{$nam}){push(@order,$nam);}
			$STO{$nam} .= $seq;
			$s++;
		}
	}
}
close(STO);
my @MA = split(//,$STO{$order[0]});
for my $nam (@order)
{
	my $seq;
	my $n = 0;
	while(exists $MA[$n])
	{
		if($MA[$n] ne "-"){$seq .= substr($STO{$nam},$n,1);}
		$n++;
	}
	print ">".$nam."\n";
	print $seq."\n";
}
