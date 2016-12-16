#!/usr/bin/perl

open(A,$ARGV[0]); # Network.csv from accnet
@txt = <A>;
close A;

shift(@txt);

foreach $l (@txt)
{
	@c = split(/\t/,$l);
	$type1{$c[0]}++;
	$type2{$c[1]}++;
	$matrix{$c[0]}{$c[1]} = 1;
}

@strains = sort(keys(%type2));
$names = join("\t",@strains);

print "$names\n";
foreach $k (keys(%type1))
{
	print "$k\t";
	foreach $j (@strains)
	{
		if (exists($matrix{$k}{$j}))
		{
			print "1\t";
		}else{
			print "0\t";
		}
	}
	print "\n";
}
