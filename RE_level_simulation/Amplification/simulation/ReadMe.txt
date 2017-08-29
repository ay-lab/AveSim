Requirement :
        R/perl environment
        R packages : splitstackshape and BSgenome.Hsapiens.UCSC.hg19

amp.simulation.config.txt file is need to run the simulation. This file contains the following the information

amp_chrom = chr2 #chromosome where amplification event will be performed
amp_start = 1000000 #Starting coordinate of amplification
amp_end   = 2000000 #Ending coordinate of amplification
amp_put   = 2000001 #Coordinate where the amplified region should be placed 
cnv_copy  = 1 #How many copies should be there. With very high number this will represent a Homogeneously staining region simulation
hetero    = 0.50 #heterogeneity level i.e. Final fastq file will contain 50% reads from amplified chromosome and 50% from normal chromosome
sim_rest  = no #Should all chromosomes be simulated (yes/no). yes will produce all the 23 chromosme fastq files.
reads     = 5000000 #Each chromosome read number

To run, type "perl simulate_amplification.pl". 
