Requirement :
        R/perl environment
        R packages : splitstackshape and BSgenome.Hsapiens.UCSC.hg19

trans.simulation.config.txt file is need to run the simulation. This file contains the following the information

don_chrom = chr2 #Donor chromosome
don_start = 1000000 #Starting position of the translocated region
don_end   = 2000000 #Ending position of the translocated region
rec_chrom = chr3 #Receiving chromosome
rec_put   = 2000000 #Coordinate to place the donor chromosome region. Note: (1-1999999) region will be deleted in the chromosome 3
hetero    = 0.50 #heterogeneity level i.e. Final fastq file will contain 50% reads from translocated chromosome and 50% from normal chromosome
sim_rest  = no #Should all chromosomes be simulated (yes/no). yes will produce all the 23 chromosme fastq files.
reads     = 5000000 #Each chromosome read number

To run, type "perl simulate_translocation.pl".
