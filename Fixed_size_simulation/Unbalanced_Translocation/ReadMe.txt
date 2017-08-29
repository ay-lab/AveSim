Requirements
R package gplots
can be installed by install.packages("gplots") command from R interactive enviornment.

To run please do
perl masterScript.pl TransFrameWork.txt

The "TransFrameWork.txt" file has all the information to generate a inter-chromosomal Hi-C matrix will many different types of translocation types. The file description is as follows,

#Input Files
bedFile ../InputFiles/IMR90.chr19_chr20.bed	#Bed file of the Hi-C matrix containing the indices
matFile ../InputFiles/IMR90.chr19_chr20.matrix #matrix file in sparse matrix format
################### 
#Output File Prefix
output  IMR90
###################
#Chromosome sizes (in 40Kb bins)
chrA	1	1479  #start 	End
chrB	1480	3055  #Start	End
###################
#1st Chromosome and Its changes in Intra-chromosomal matrix (chrA)
chr	chr19 #1st chromosome name
###################
#2nd Chromosome and Its changes Intra-chromosomal matrix (chrB)
chr	chr20 #Second chromosome name
###################
#Contact count treatment
chrA	HMD #After translocation, how do you want to treat the chromosome A intra-chromosomal contacts: HMD/Homogeneous or HTD/Hetergeneous deletion
chrB	HMD #Same as above but in Chromosome A 
###################
#Translocation regions in Inter-chromosomal matrix
This section controls how and what type translocations is required. Different types can be described like

For reciprocal balanced translocation like chrA + chrB = chrA1~chrB2) + chrB1~chrA2)
e.g.
New_Chromosome1	1:591~1953:3055)
New_Chromosome2	1480:1952~592:1479)

For Non-reciprocal balanced translocation like chrA + chrB = chrA1) + chrB1~chrA2)chrB2
e.g.
New_Chromosome1 1:591)
New_Chromosome2 1480:2110~592:1479)2111:3055)

For Unbalanced translocation like chrA + chrB = chrA1) + chrB1~chrA2)
e.g.
New_Chromosome1 1:591)
New_Chromosome2 1480:2110~592:1479)
###################
#Contact counts
Intra	self	#Keep the original intra-chromosomal count as it is
Inter	self #Use intra-chromosomal contacts to calculate mean inter-chromosomal counts. 
################### 
