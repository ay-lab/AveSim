library(BSgenome.Hsapiens.UCSC.hg19)
args <- commandArgs(trailingOnly = TRUE) 
#Rscript deletion.r chr1 150000000 240000000
chr = read.table("../scripts/chr.list",header=F)
if (as.character(args[4]) == "yes"){
	j = 1
	while (j <= length(chr$V1)){
		if (as.character(args[1]) != as.character(chr$V1[j])){
			x.a = getSeq(Hsapiens, as.character(chr$V1[j]))
			sequence = as.character(x.a)
			name = paste(">",chr$V1[j],sep="")
			fileConn <- file(as.character(paste("Normal.",chr$V1[j],".genome.fasta",sep="")))
			header = as.character(paste(name,sequence,sep="\n"))
		        writeLines(header,fileConn)
			close(fileConn)
		}
		j = j+1
	}
}
full.A = getSeq(Hsapiens, as.character(args[1]))
a = getSeq(Hsapiens, as.character(args[1]),start=start(full.A)[[1]],end=as.numeric(args[2])-1)
b = getSeq(Hsapiens, as.character(args[1]),start=as.numeric(args[2]), end=as.numeric(args[3]))
c = getSeq(Hsapiens, as.character(args[1]),start=as.numeric(args[3])+1, end=end(full.A)[[1]])
a = as.character(a)
b = as.character(b)
c = as.character(c)
fileConn <- file(as.character(paste("Normal.",args[1],".genome.fasta",sep="")))
sequence = as.character(paste(full.A,sep=""))
name = paste(">",args[1],sep="")
header = as.character(paste(name,sequence,sep="\n"))
writeLines(header,fileConn)
close(fileConn)

fileConn <- file(as.character(paste("Deleted.",args[1],".genome.fasta",sep="")))
sequence = as.character(paste(a,c,sep=""))
name = paste(">",args[1],sep="")
header = as.character(paste(name,sequence,sep="\n"))
writeLines(header,fileConn)
close(fileConn)
