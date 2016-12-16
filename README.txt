AccNET V1.1.3 - Accessory Constellation Network.

Last update: 09/03/2016

Developed by: Val F. Lanza. (valfernandez.vf@gmail.com)

DESCRIPTION

AccNET is a comparative genomic tool for accessory genome analysis using
bipartite networks. The software has been designed to be compatible with
most of the Network Analysis software (i.e. Cytoscape, Gephi or R).

AccNET has been developed in Perl and it is designed for Linux 
platforms. Please read the Dependencies secction for more details. The 
software builds a bipartite network integrated by two kind of nodes 
"Genomic Units (GU)" and "Homologous Proteins Cluster (HPC)". GU can be single 
elements such chromosomes or plasmids, or complex set such as genomes, 
pangenomes or even enviromental proteomes.


INPUT DATA

AccNET works with proteomes. Each proteome must be in a single file.
AccNET do not works with DNA data. A proteome can be a single element 
such as Chromosome, plasmid, phage etc... or complex element (Genome 
with a mix of chromosome and plasmids) but in any case, each element is 
defined by its file.

OUTPUT DATA

 -Network.csv:	      This is the network definition and include three 
				      columns: "Source", "Target", "Weigth" and "Type".
 -Table.csv:	      This file include all nodes attribute information.
	
 -Representatives.faa FASTA file with representative AA sequence of
					  each cluster (HPC).

	please read the VISUALIZATION secction.



EXAMPLES

Accesory Network for genomes.

	Simple:		accnet.pl --in *.faa
	Advance:	accnet.pl --in *.faa --threshold 0.8 
						  --kp '-s 1.5 -e 1e-8 -c 0.8' 
						  --out Network_example.csv 
						  --tblout Table_example.csv
						  --fast yes
	
Whole genomes. Only recommended for plasmids or inter-species comparisson.

				accnet.pl --in *.faa --threshold 1.1

VISUALIZATION:

#Gephi visualization (https://gephi.org/).
	-Open Gephi.
	-Make a new Project. (File -> New Project)
	-Import spreadsheet (File -> Import spreadsheet...)
	-Select "Network.csv" as "Edges Table"
	-Import spreadsheet (File -> Import spreadsheet...)
	-Select "Table.csv" as "Nodes Table"

#Cytoscape visualization (http://www.cytoscape.org/)
	-version 2.8.x
		-Import Network file (File -> Import -> Network from Table)
		-Select "Network.csv"
		-Remove 1st line ("Show Text File Import Options" 
							-> "Transfer first line as attribute names")
		-Select delimiter "Tab"
		-Select 1st column as "Source Interaction"
		-Select 2nd column as "Target Interaction"
		-Check "Weight" column to import.
		-Import.
		
		-Import Node Attributes (File -> Import -> Attibutes from Table)
		-Select "Table.csv" file
		-Select delimiter "Tab"
		-Import column headers ("Show Text File Import Options"
							-> "Transfer first line as attribute names")
		-Import
	
	-version 3.x
		-Import Network file (File -> Import -> Network -> File)
		-Select "Network.csv"
		-Remove 1st line ("Show Text File Import Options"
							->"Transfer first line as attribute names ")
		-Select delimiter "Tab"
		-Select 1st column as "Source Interaction"
		-Select 2nd column as "Target Interaction"
		-Check "Weight" column to import.
		-Import.
		
		-Import Node Attributes (File -> Import -> Table -> File)
		-Select "Table.csv" file
		-Select delimiter "Tab"
		-Import column headers ("Show Text File Import Options"
							 ->"Transfer first line as attribute names")
		-Import
			

NETWORK CLUSTERING

Since AccNET v1.1.3 Clustering.r script has been added to the project. 
Clustering.r script performs a clustering analysis that found both GU 
and HpC clusters based on the network adjacent matrix. Clustering.r are 
written in R language and requires the libraries dplyr, tidyverse, cluster
and mclust. GU clusters are calculated by two methods: first with mclust 
(Gaussian Mixture Modelling for Model-Based Clustering,Classification, and 
Density Estimation) and second by hierarchical clustering. In HpC case, 
the clusters are only calculated from hierarchical clustering method. 
Both methods, hierarchical and bayesian use a distance matrix as input data.
This distance matrix are calculated using the distance binary method. In GU 
case, the GU are taken as objects and HpC as variables and vice versa in HpC
case. For hierarchical clustering different heights are taken to create the
clusters. The cut points are calculated as the quantiles 75, 85, 90, 95 and
99 of tree heights. The resulting output file is a tab format file that 
can be loaded in Gephi or Cytoscape.

Command line example (command has to be executed in the Network.csv folder):

	R -f /path/to/Clusering.r

Installing dependencies:
Open R and type:
	install.packages(tidyverse)
	install.packages(cluster)
	install.packages(mclust)
			
			
			
DEPENDENCIES;

All dependencies must be included in PATH variable.

	- kClust: ftp://toolkit.lmb.uni-muenchen.de/pub/kClust/
	  (Hauser et al., 2013) kClust: fast and sensitive clustering 
	  of large protein sequence databases. BMC bioinformatics.
	  PMID: 23945046

	- trimAl: http://trimal.cgenomics.org/downloads
	  trimAl: a tool for automated alignment trimming in large-scale 
	  phylogenetic analyses.Salvador Capella-Gutierrez; Jose M. Silla-Martinez; 
	  Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.
	
	- muscle: http://www.drive5.com/muscle/ & linux common repositories
	  Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high 
	  accuracy and high throughput. Nucleic Acids Res. 32(5):1792-1797
	  
	  *Linux command line installation: "sudo apt-get install muscle" or 
							"sudo yum install muscle"

	- PHYLIP: http://evolution.genetics.washington.edu/phylip.html 
			  & linux common repositories
	  Felsenstein, J. (1993). PHYLIP (phylogeny inference package), 3.5 c ed. 
	  Department of Genetics, University of Washington.
	  
	  *Linux command line installation: "sudo apt-get install phylip" or 
							"sudo yum install phylip"
	  
	Perl packages dependencies:
		-List::Util  		(Core-modules)
		-Getopt::Long 		(Core-modules)
