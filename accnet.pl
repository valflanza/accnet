#!/usr/bin/perl

use List::Util; qw(first max maxstr min minstr reduce shuffle sum);
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;
use warnings;
use Statistics::R;


  
  




#### Variables

my $indexStrain;
my $indexProt;
my $tax;
my $output;
my $finalOutput;
my $outFile;
my $outTable;
my $kClustParameters = "-s 1.12 -c 0.8 -e 1e-4";  ###options by default;
my $progress;
my $clean;
my $threshold;
my $fast;
my $total;
my $clustering;



my @text;
my @allFasta;
my @cls;
my @headers;
my @net;
my @salida;
my @inFiles;
my @Options;
my @inNet;

my %taxasHash;
my %referenceHash;
my %sequences;
my %fasta;
my %indexHeader;
my %clusters;
my %clustersMembers;
my %Synt;


###temp variables
my @c;
my @c2;
my $l;
my $ar;
my $k;
my $j;
my $i;
my $n;
my $tmp;
my $actual;
my $firstLine;
my $map;
my $file;
my $strain;
my $alnIn;
my $alnOut;
my %suma;
my @dists;
my $dist;
my @txt;
my @val;
my $aln;
my $member;
my @uniqs;
my %repeat;
my @cSynt;
my %clusterTwin;
my @tmpArray;
my %hashTwinGroup;
my %maxCls;
my $tmpHeader;
my $tmpCluster;
my $id;
my @rename;
my %renameHash;
my @clusterFasta;
my $actualNumCluster;
my $R;

my $PATH = abs_path($0);
$PATH =~ s/\/accnet.pl//;
my $RUNPATH = getcwd();


############### SET COMMAND LINE OPTIONS #################################

@Options = (
		
		{OPT=>"in=s{,}",	VAR=>\@inFiles,	DESC=>"Proteome Files to analyze"},
		{OPT=>"out=s",	VAR=>\$outFile,	DEFAULT => 'Network.csv', DESC=>"Network filename"},
		{OPT=>"kp=s",	VAR=>\$kClustParameters,	DEFAULT => '-s 1.12 -c 0.8 -e 1e-4', DESC=>"'kClust parameters'"},
		{OPT=>"tblout=s",	VAR=>\$outTable,	DEFAULT => 'Table.csv', DESC=>"Table filename"},
		{OPT=>"threshold=s",	VAR=>\$threshold,	DEFAULT => '1', DESC=>"Percent of genomes to consider coregenome (values > 1 includes coregenome to network)"},
		{OPT=>"fast=s",	VAR=>\$fast,	DEFAULT => 'no', DESC=>"Skip the phylogenetic distance determination"},
		{OPT=>"clustering=s", VAR=>\$clustering, DEFAULT => 'no', DESC=>"Perform the network clustering process?"},
		{OPT=>"clean=s",	VAR=>\$clean,	DEFAULT => 'yes', DESC=>"Remove all the intermediary files"}
		
	);


#Check options and set variables
(@ARGV < 1) && (usage());
GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

# Now setup default values.
foreach (@Options) {
	if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
	${$_->{VAR}} = $_->{DEFAULT};
	}
}
	

unless(@inFiles){
	print STDERR "You must specified the input files\n";
	&usage();
}
	



###########################  InputFile renaming
$indexStrain=1;
open(O,">equivalence.txt");
foreach $ar (@inFiles)
{
	@c2 = split("/",$ar);
	@c = split('\.',$c2[-1]);
	
	$tax=$c[0];
	$tax =~ s/$RUNPATH//;
	$taxasHash{$indexStrain}=$tax;
	
	$indexProt=0;
	
	
	open(A,$ar);
	@text = <A>;
	close A;
	
	foreach $l (@text)
	{
		if($l =~ /^>/)
		{
			$indexProt++;
			$referenceHash{$indexStrain}{$indexProt} = $l; ### Save protein information
			$referenceHash{$indexStrain}{$indexProt} =~ s/>//;
			$sequences{$indexStrain}{$indexProt} = ">$indexStrain|$indexProt\n";  ##New header
			print O ">$indexStrain|$indexProt\t$l";
		}else{
			$sequences{$indexStrain}{$indexProt} .= $l;
		}
	}
	$indexStrain++;
}
close O;

open(O,">all_fasta.tmp");

foreach $k (keys(%sequences))
{
	foreach $j (keys(%{$sequences{$k}}))
	{
		print O "$sequences{$k}{$j}";
	}
}
close O;
######################## End InputFile renaming

########## Clustering and post-proccess
system("$PATH/bin/kClust -i all_fasta.tmp -d ./kclust $kClustParameters");

open(FASTA,"all_fasta.tmp");
@allFasta = <FASTA>;
close FASTA;

foreach $l (@allFasta) ## hashing fasta file ##
{
	if($l =~ />/)
	{
		$actual = $l;
		$actual =~ s/\s//g;
	}
	$fasta{$actual}.= $l;
}



open(HEAD,"./kclust/headers.dmp");
@headers = <HEAD>;
close HEAD;

##### headers.dmp ###############
#	1 >NC_014534|28		#
#	88 >NC_009930|107	# File headers.dmp (index_sequence -> header)
#	5301 >NC_010631|164	#
#################################

foreach $l (@headers)
{
	chomp $l;
	@c = split(' ',$l);
	$c[1] =~ s/\s//g;
	$indexHeader{$c[0]} = $c[1];
	
}

open(CLS,"./kclust/clusters.dmp");
@cls = <CLS>;
close CLS;

shift(@cls);
############ clusters.dmp ###############
#	# 40086 ->header (must remove)	#
#	1 1				#
#	88 1				#
#	5301 1				# File clusters.dmp (index_sequence -> index_cluster(same as representative_sequence))
#	20249 1				#
#	40027 1				#
#	40068 1				#
#	2 2				#
#########################################

foreach $l (@cls)
{
	chomp $l;
	@c = split(' ',$l);
	push(@{$clusters{$c[1]}},$fasta{$indexHeader{$c[0]}});
	push(@{$clustersMembers{$c[1]}},$indexHeader{$c[0]});
	
}
################### END Clustering and post-proccess





################## Distances calculation


$actualNumCluster=0;
$total = scalar(keys(%clusters));
$| = 1;  # Turn off buffering on STDOUT.

foreach $k (keys(%clusters))
{
	
	$output = join('',@{$clusters{$k}});
	open(OUT,">Cluster_$k.fasta");
	print OUT $output;
	close OUT;
	undef @uniqs;
	undef %repeat;
	foreach $member (@{$clustersMembers{$k}}) 
	{
		if(exists($repeat{$member}))
		{
			$repeat{$member}++;
		}else{
			push(@uniqs,$member);
			$repeat{$member}++;
		}
	}
	
	if (scalar(@uniqs) < scalar (@inFiles) * $threshold)
	{
		push(@inNet,$k);
		if(scalar(@{$clusters{$k}}) <2)
		{
			open(TMP,"Cluster_$k.fasta");
			$firstLine = <TMP>;
			close TMP;
			$firstLine =~ s/>//;
			@c = split('\|',$firstLine);
			push(@net,"Cluster_$k\t$taxasHash{$c[0]}\t".scalar(@inFiles)/2 . "\tUndirected");
		}else{
			if($fast eq 'Y' | $fast eq 'Yes' | $fast eq 'yes' | $fast eq 'y')
			{
				open(TMP,"Cluster_$k.fasta");
				@tmpArray = <TMP>;
				@clusterFasta = grep(/>/,@tmpArray);
				foreach $l (@clusterFasta)
				{
					$l =~ s/>//;
					@c = split('\|',$l);
					push(@net,"Cluster_$k\t$taxasHash{$c[0]}\t".scalar(@inFiles)/2 . "\tUndirected");
				}
				
				close TMP;
			}else{
				push(@net,distance("Cluster_$k"));
			}
		}
	}
	$actualNumCluster++;
	
	print("\rProccessing cluster $actualNumCluster of $total");
	

}
print "\n";

################# Output section



###Output Network
$finalOutput = join("\n",@net);
open(FINAL,">$outFile");
print FINAL "Source\tTarget\tWeight\tType\n";
print FINAL $finalOutput;
close FINAL;


###Output Annotation Table

%Synt = synteny();

open(TABLE, ">$outTable");
print TABLE "ID\tType\tTwinGroup\tDescription\n";
foreach $k (keys(%taxasHash))
{
	print TABLE "$taxasHash{$k}\tGU\tGU\t$k\n";
}

foreach $k (@inNet)
{
	$tmp = $indexHeader{$k};
	$tmp =~ s/>//;
	@c = split('\|',$tmp);
	
	$id = "Cluster_$k";
	if (exists($Synt{$id}))
	{
		print TABLE "$id\tCluster\t$Synt{$id}\t$referenceHash{$c[0]}{$c[1]}";
	}else{
		print TABLE "$id\tCluster\tSingle\t$referenceHash{$c[0]}{$c[1]}";
	}
		 
}
close TABLE;

#output representative.fasta

system("grep '>' Cluster*.fasta | sed 's/.fasta:/\t/' > rename.tmp");
open(A,"rename.tmp");
@rename = <A>;
close A;

foreach $l (@rename)
{
	chomp $l;
	@c = split('\t',$l);
	$renameHash{$c[1]} = $c[0];
}
open(FASTA,"./kclust/representatives.fas");
open(OUT, ">Representatives.faa");
while ($l = <FASTA>)
{
	if ($l =~ />/)
	{
		chomp $l;
		print OUT ">$renameHash{$l}\n";
	}else{
		print OUT $l;
	}
}
close FASTA;


if($clustering eq 'Y' | $clustering eq 'Yes' | $clustering eq 'yes' | $clustering eq 'y')
{
	clusteringNetwork();
}

if($clean eq "yes")
{
	system("rm -r Cluster_* kclust outfile seq_temp temp all_fasta.tmp");
}







######################## Subrutines #############################
sub synteny
{
	
	
	foreach $l (@net)
	{
		chomp $l;
		@c = split('\t',$l);
		push(@{$clusterTwin{$c[0]}},$c[1]);
	
	}
	foreach $k (keys(%clusterTwin))
	{
		if (@{$clusterTwin{$k}}>1)
		{
			@tmpArray = sort(@{$clusterTwin{$k}});
			$tmpHeader = join('=',@tmpArray);
			push(@{$maxCls{$tmpHeader}},$k);
		}
	}

	$i=0;
	foreach $k (keys(%maxCls))
	{
		if (scalar(@{$maxCls{$k}}>1))
		{
			$i++;
			foreach $tmpCluster (@{$maxCls{$k}})
			{
				
				$hashTwinGroup{$tmpCluster} = "Twin_$i";
			
			}
		}
	}
	return %hashTwinGroup;
}

sub distance 
{
	if( -e "outfile")
	{
		system("rm outfile");
	}
	if( -e "seq_temp")
	{
		system("rm seq_temp");
	}
	if ( -e "temp")
	{
		system("rm temp");
	}
	
	$file = $_[0];
	
	undef @salida;
	undef $map;
	undef %suma;
	undef @dists;
	undef $alnIn;
	undef $alnOut;
	undef $aln;
	undef @txt;
	
	
	
	
	
	system("$PATH/bin/muscle -quiet -in $file.fasta -out $file.aln");
	
	system("$PATH/bin/trimal -strictplus -phylip -in $file.aln -out $file.aln.phylip >/dev/null 2>&1");
###### Protdist Process 


	open(TMP,">seq_temp");
	print TMP "$file.aln.phylip\nY\n";
	close TMP;
	system("$PATH/bin/protdist < seq_temp > temp");
	
	#system("cat outfile");
	if(open(D,"outfile"))
	{@txt = <D>;}
	close D;

	
	shift(@txt); ##remove header	

	for($i=0; $i<scalar(@txt); $i++)
	{
		if($txt[$i] =~ /^\d/)
		{
			@c = split(' ',$txt[$i]);
			@c2 = split('\|',$c[0]);			
			$actual = $c2[0];
			$suma{$actual} =0;
	
			for($j=1;$j<scalar(@c);$j++)
			{
				$suma{$actual} += abs($c[$j]);
			}
			
		}else{
			@c = split(' ',$txt[$i]);
			for($j=0;$j<scalar(@c);$j++)
			{
				$suma{$actual} += abs($c[$j]);
			}
		}
	}
	@val = values(%suma);
	@dists =  sort {$a <=> $b} @val;
	
	foreach $strain (keys(%suma))
	{
		$n = scalar(keys(%suma));
		if($suma{$strain} eq 0)
		{
			$dist = $n/2;
		}else{
			$dist = log(($n*$n)/$suma{$strain});
			if($dist <= 0)
			{
				$dist =0.01;
			}
		}
		
		push(@salida,"$file\t$taxasHash{$strain}\t$dist\tUndirected");  

	}
	
	return @salida;
}	
######Clustering Network subrutine

sub clusteringNetwork {
	$R = Statistics::R->new(shared => 1);
	$R->start();
	$R->set('TableFile',$outFile);
	$R->set('runPath',$RUNPATH);
	$R->run_from_file("$PATH/lib/Clustering.r");
	$R-> stop();
}

######Usage subrutine
sub usage {
	print "Usage:\n\n";
	print "Accesory genes Network (ACCNET) v1.2\n";
	print "writen by: Val F. Lanza (valfernandez.vf\@gmail.com)\n\n";
	print "Accesory Network for genomes.\n\nSimple:		accnet.pl --in *.faa\nAdvance:	accnet.pl --in *.faa --threshold 0.8 --kp '-s 1.5 -e 1e-8 -c 0.8' --out Network_example.csv --tblout Table_example.csv --fast Yes --clustering Yes\n\nWhole genomes. Only recommended for plasmids or inter-species comparisson.\n\n		accnet.pl --in *.faa --threshold 1.1\n";

	print "\nParameters:\n\n";
	foreach (@Options) {
		
		printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
			defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
	}
	print "\n\n\n";
	exit(1);
}

