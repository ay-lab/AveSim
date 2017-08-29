Requirement :
	R/perl environment
	R packages : splitstackshape and BSgenome.Hsapiens.UCSC.hg19

del.simulation.config.txt file is need to run the simulation. This file contains the following the information

del_chrom = chr2 #chromosome where deletion event will be performed
del_start = 1000000 #Starting coordinate of deletion
del_end   = 2000000 #Ending coordinate of deletion
hetero    = 0.50 #heterogeneity level i.e. Final fastq file will contain 50% reads from deleted chromosome and 50% from normal chromosome.
sim_rest  = no #Should all chromosomes be simulated (yes/no). yes will produce all the 23 chromosme fastq files.
reads     = 5000000 # Each chromosome read number 

To run, type "perl simulate_deletion.pl".
