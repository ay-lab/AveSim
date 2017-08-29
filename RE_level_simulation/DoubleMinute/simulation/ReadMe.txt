Requirement :
        R/perl environment
        R packages : splitstackshape and BSgenome.Hsapiens.UCSC.hg19

dm.simulation.config.txt file is need to run the simulation. This file contains the following the information

dm_chrom = chr2 #chromosome where double minute event will be performed
dm_start = 2000000 #Starting coordinate of dm
dm_end   = 5000000 #Ending coordinate of dm
reads    = 5000000 #read read number

To run, type "perl simulate_doubleminute.pl". This will create Normal_DM.re_R1.fastq and Normal_DM.re_R2.fastq as the final fastq files. 
