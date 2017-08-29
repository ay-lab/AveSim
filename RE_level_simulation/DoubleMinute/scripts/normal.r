library(BSgenome.Hsapiens.UCSC.hg19)
args <- commandArgs(trailingOnly = TRUE) 
#Rscript duplication chr1 chr2 150000000 240000000 3000000 3000001 1
chr = read.table("../scripts/chr.list",header=F)
j = 1
while (j <= length(chr$V1)){
	x.a = getSeq(Hsapiens, as.character(chr$V1[j]))
	sequence = as.character(x.a)
	name = paste(">",chr$V1[j],sep="")
	fileConn <- file(as.character(paste("Normal.",chr$V1[j],".genome.fasta",sep="")))
	header = as.character(paste(name,sequence,sep="\n"))
        writeLines(header,fileConn)
	close(fileConn)
	j = j+1
}
