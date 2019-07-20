#!/usr/bin/perl
########################################################################################################################
# Stitching alignments v1.6
# v1.1 - adding support for uni2loc_fas: this a pre-annotated uniprot sequence database that can be used with jackhmmer
#        http://robetta.bakerlab.org/downloads/contact_prediction/uni2loc_2014_03_18.fas.zip
# v1.2 - fixed bug to deal with alignments that contain multiple hits with same header
# v1.3 - bug fix related to v1.2 fix
# v1.4 - fixed a bug when dealing with the new uniport id(s), for jackhmmer pre-annotated database
# v1.5 - fixed a bug when dealing with the new uniport id(s), for hhblits
# v1.6 - fixed a major bug when dealing with the new uniport id(s) that are 10 characters in length.
#      - Thanks to Abantika Pal for catching this!
########################################################################################################################

# Email: krypton@uw.edu (please contact me to report bugs or issues)
# http://gremlin.bakerlab.org/complexes

# If you find the following script/method useful, please cite us in your work
# Sergey Ovchinnikov, Hetunandan Kamisetty, David Baker
# Robust and accurate prediction of residue-residue interactions across protein interfaces using evolutionary information
# Elife (2014)
# DOI: http://dx.doi.org/10.7554/eLife.02030

# The following paths need to be defined if you are going to use "-uni2loc" option
my $sgrep = "/work/krypton/bin/sgrep"; # sgrep (sorted grep): http://sgrep.sourceforge.net/

# Latest version can be obtained here: http://gremlin.bakerlab.org/uni2loc/
my $uni2loc_db = "/scratch/USERS/database/uni2loc_2014_03_18";

unless(exists $ARGV[0])
{
	print "=================================================================\n";
	print "           Stitching alignments v1.6   \n";
	print "=================================================================\n";
	print " -a            a3m/fas of alignment A [REQUIRED]\n";
	print " -b            a3m/fas of alignment B [REQUIRED]\n";
	print " -min          delta_gene [default = 01]\n";
	print " -max          delta_gene [default = 20]\n";
	print " -out          location for fasta out file [OPTIONAL]\n";
	print "-----------------------------------------------------------------\n";
	print " Use following tags to enable options\n";
	print "-----------------------------------------------------------------\n";
	print " -strand       only stitch genes from same strand-ness F-F vs R-R\n";
	print "               one of the uni2loc options are required\n";
	print " -uni2loc      use uni2loc database (specified inside script)\n";
	print " -uni2loc_fas  for alignments generated using uni2loc_fas db\n";
	print "               otherwise a3m input from from hhblits assumed\n";
	print "=================================================================\n";
	die();
}

my @a3m;
while ($arg = shift())
{
	if ($arg =~ m/^-/)
	{
		$arg = lc($arg);
		if ($arg eq "-a")		{$a3m[0] = shift(); next;}
		if ($arg eq "-b")		{$a3m[1] = shift(); next;}
		if ($arg eq "-min")		{$min = shift(); next;}
		if ($arg eq "-max")		{$max = shift(); next;}
		if ($arg eq "-out")		{$out = shift(); next;}
		if ($arg eq "-uni2loc")		{$uni2loc = 1; next;}
		if ($arg eq "-uni2loc_fas")	{$uni2loc_fas = 1; next;}
		if ($arg eq "-strand")		{$strand = 1; next;}
	}
}
if(defined $uni2loc){for my $n ($sgrep,$uni2loc_db){unless(-e $n){die("'$n' not found");}}}

unless(-e $a3m[0]){die("-a '$a3m[0]'");}
unless(-e $a3m[1]){die("-b '$a3m[1]'");}
unless(defined $min){$min = 1;}
unless(defined $max){$max = 20;}

##########################
# convert the uniprot accession_id into a number
# http://www.uniprot.org/manual/accession_numbers
# [5]       [4]   [3]       [2]       [1]       [0]   [3]   [2]       [1]       [0]
# [A-N,R-Z] [0-9] [A-Z]     [A-Z,0-9] [A-Z,0-9] [0-9] [A-Z] [A-Z,0-9] [A-Z,0-9] [0-9]
# [A-N,R-Z] [0-9] [A-Z]     [A-Z,0-9] [A-Z,0-9] [0-9]
# [O,P,Q]   [0-9] [A-Z,0-9] [A-Z,0-9] [A-Z,0-9] [0-9]

my (%pa,%ma); for $m (A..N,R..Z){$pa{$m} = "A";}for $m (O..Q){$pa{$m} = "B";}
my $n = 0;for my $t (0..9){$ma{A}[0]{$t} = $n; $ma{A}[4]{$t} = $n;$ma{B}[0]{$t} = $n; $ma{B}[4]{$t} = $n;$n++;}
my $n = 0;for my $t (A..Z,0..9){$ma{A}[1]{$t} = $n; $ma{A}[2]{$t} = $n;$ma{B}[1]{$t} = $n; $ma{B}[2]{$t} = $n; $ma{B}[3]{$t} = $n;$n++;}
my $n = 0;for my $t (A..Z){$ma{A}[3]{$t} = $n;$ma{A}[5]{$t} = $n;$ma{B}[5]{$t} = $n;$n++;}
sub uni_num
{
	my $uni = substr($_[0],0,6);
	my $p = $pa{substr($uni,0,1)};
	my ($tot,$num) = (1,0);
	if(length($_[0]) == 10) # for the extended uniprot accession see: http://www.uniprot.org/changes/accession_format
	{
		my $u = reverse(substr($_[0],6,4)); #bug fix v1.6
		for my $n (0..3)
		{
			$num += $ma{$p}[$n]{substr($u,$n,1)}*$tot;
			$tot *= int(keys %{$ma{$p}[$n]});
		}
	}
	my $u = reverse($uni);
	for my $n (0..5)
	{
		$num += $ma{$p}[$n]{substr($u,$n,1)}*$tot;
		$tot *= int(keys %{$ma{$p}[$n]});
	}
	return $num;
}
my @nums = (0..9,a..z,A..Z);
my %nums = map { $nums[$_] => $_ } 0..$#nums;
sub base2num
{
	my ($base,$rep) = @_;my $number = 0;
	for( $rep =~ /./g ){$number *= $base;$number += $nums{$_};}
	return $number;
}
my %seq;
my %SEEN;
my $d = 0;
my @DO;
for my $i (0,1)
{
	my $add_seq = 1;
	my $k = 0;
	my $uni = "null";
	open(A3M,$a3m[$i]);
	while($line = <A3M>)
	{
		chomp($line);
		if(substr($line,0,1) eq ">")
		{
			if($k > 0) # not the first sequence
			{
				if(defined $uni2loc_fas)
				{
					# assuming alignment is in fasta format, where the header looks as follows:
					# >H4SB28_3o4Y_001d_F/18-478
					my @tm = split(/[\_\/\s]/,substr($line,1)); #bug fix v1.4
					$uni = $tm[0];
					@{$DO[$d]} = (base2num(62,$tm[2]),$tm[1],$i,$uni.$i,$tm[3]);$d++;
				}
				else
				{
					# assuming alignment is coming from HHsuite, where the header looks as follows:
					# >tr|E8G9Y7|E8G9Y7_SALMO NADH dehydrogenase subunit J OX=882866 OS=IA_2010008284. GN=SEEM8284_18765 PE=3 SV=1
					# >tr|A0A0A0MTS7|A0A0A0MTS7_HUMAN Titin OS=Homo sapiens GN=TTN PE=1 SV=1

					my @arr = split(/[\|\s]+/,$line);
					if(exists $arr[1] and (length($arr[1]) == 6 or length($arr[1]) == 10)) #bug fix v1.5
					{
					
						my %info; for my $tm (@arr){my @ar = split("=",$tm);if(exists $ar[1]){$info{$ar[0]} = $ar[1];}}
						$uni = $arr[1];	

						if(defined $uni2loc)
						{
							for my $loc (`$sgrep $uni $uni2loc_db`)
							{
								chomp($loc);my @tm = split(/[\t\_]+/,$loc);
								if($tm[0] eq $uni) #bug fix v1.4
								{
									@{$DO[$d]} = (base2num(62,$tm[2]),$tm[1],$i,$uni.$i,$tm[3]);$d++;
								}
							}
						}
						# if no uni2loc database is provided "OX" information is used to determine if the 
						# two genes are from the same genome.
						elsif(defined $info{OX}){@{$DO[$d]} = (uni_num($uni),$info{OX},$i,$uni.$i,"X"); $d++;}
						else{print STDERR "NO_OX_DEFINITION: $line\n";}
					}
					else{print STDERR "NOT_VALID_UNIPROT_HHBLITS_HEADER: $line\n";}
				}
				#bug fix v1.2: for duplicates
				unless(exists $SEEN{$uni.$i}){$SEEN{$uni.$i} = 1; $add_seq = 1;}else{$add_seq = 0;}
			}
		}
		elsif($add_seq == 1)
		{
			$line =~ s/[a-z]//g;
			$seq{$uni.$i} .= $line;
		}
		$k++;
		
	}
	close(A3M);
}
@DO = sort { $a->[1] cmp $b->[1] } @DO;
my %NET;
my $i = 0;
while(exists $DO[$i])
{
	my $j = $i + 1;
	while(exists $DO[$j] and $DO[$i][1] eq $DO[$j][1]) #if same genome
	{
		if($DO[$i][2] ne $DO[$j][2]) #if different alignments
		{
			my $diff = abs($DO[$i][0]-$DO[$j][0]); #delta_gene

			#start building network 
			if((!exists $NET{$DO[$i][3]} or $diff < $NET{$DO[$i][3]}[0]) and (!defined $strand or $DO[$i][4] eq $DO[$j][4])){@{$NET{$DO[$i][3]}} = ($diff,$DO[$j][3]);}
			if((!exists $NET{$DO[$j][3]} or $diff < $NET{$DO[$j][3]}[0]) and (!defined $strand or $DO[$i][4] eq $DO[$j][4])){@{$NET{$DO[$j][3]}} = ($diff,$DO[$i][3]);}
		}
		$j++;
	}
	$i++;
}
if(defined $out) # OPEN FILE
{
	my ($uni_a,$uni_b) = ("null0","null1");
	open(OUT,">$out");
	print OUT ">".$uni_a."_".$uni_b."\n".$seq{$uni_a}.$seq{$uni_b}."\n";
}
my %g;
my $g_inf = 0;
my @count;
foreach my $uni_a ( keys %NET )
{
	my $uni_b = $NET{$uni_a}[1];
	if($NET{$uni_b}[1] eq $uni_a)
	{
		#keep count of $delta_gene <= 20
		if($NET{$uni_a}[0] <= 20){$count[$NET{$uni_a}[0]] += 1;} 

		if(defined $out and $NET{$uni_a}[0] >= $min and $NET{$uni_a}[0] <= $max) # WRITE TO FILE
		{
			if(substr($uni_b,-1) > substr($uni_a,-1))
			{
				print OUT ">".substr($uni_a,0,-1)."_".substr($uni_b,0,-1)."_".$NET{$uni_a}[0]."\n".$seq{$uni_a}.$seq{$uni_b}."\n";
			} 
			else
			{
				print OUT ">".substr($uni_b,0,-1)."_".substr($uni_a,0,-1)."_".$NET{$uni_a}[0]."\n".$seq{$uni_b}.$seq{$uni_a}."\n";
			}
		}
		$g_inf++;
	}
}
close(FASTA);
if(defined $out){close(OUT);} # CLOSE FILE

my $n = 0;
print "==== Distribution of delta_gene ====\n";
print "tot\t$g_inf\n";
if($g_inf > 0)
{
	my $n = 0;
	while($n <= 20)
	{
		print $n."\t".int($count[$n])."\t".sprintf("%.2f",int($count[$n])/$g_inf)."\n";
		$n++;
	}
}
