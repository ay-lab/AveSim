library(BSgenome.Hsapiens.UCSC.hg19)
args <- commandArgs(trailingOnly = TRUE) 
#Rscript duplication chr1 chr2 150000000 240000000 3000000 3000001 1
chr = read.table("../scripts/chr.list",header=F)
if (as.character(args[8]) == "yes"){
	j = 1
	while (j <= length(chr$V1)){
		if (as.character(args[1]) != as.character(chr$V1[j]) & as.character(args[2]) != as.character(chr$V1[j])){
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
if (as.numeric(args[6]) >= end(full.A)[[1]]){
        a = getSeq(Hsapiens, as.character(args[1]),start=start(full.A)[[1]],end=as.numeric(args[3])-1)
        b = getSeq(Hsapiens, as.character(args[1]),start=as.numeric(args[3]), end=as.numeric(args[4]))
        c = getSeq(Hsapiens, as.character(args[1]),start=as.numeric(args[4])+1, end=end(full.A)[[1]])
        a = as.character(a)
        b = as.character(b)
        c = as.character(c)
	b.new = paste(rep(b,as.numeric(args[7])),collapse="")
	fileConn <- file(as.character(paste("Normal.",args[1],".genome.fasta",sep="")))
        sequence = as.character(paste(full.A,sep=""))
      	name = paste(">",args[1],sep="")
       	header = as.character(paste(name,sequence,sep="\n"))
       	writeLines(header,fileConn)
       	close(fileConn)	
	fileConn <- file(as.character(paste("Duplicated.",args[1],".genome.fasta",sep="")))
        sequence = as.character(paste(full.A,b.new,sep=""))
        name = paste(">",args[1],sep="")
        header = as.character(paste(name,sequence,sep="\n"))
        writeLines(header,fileConn)
       	close(fileConn)
	}else if (as.numeric(args[6]) < end(full.A)[[1]] & as.numeric(args[5]) > start(full.A)[[1]]){
	fileConn <- file(as.character(paste("Normal.",args[1],".genome.fasta",sep="")))
	sequence = as.character(paste(full.A,sep=""))
     	name = paste(">",args[1],sep="")
       	header = as.character(paste(name,sequence,sep="\n"))
       	writeLines(header,fileConn)
       	close(fileConn)
	a = getSeq(Hsapiens, as.character(args[1]),start=start(full.A)[[1]],end=as.numeric(args[5]))
	b = getSeq(Hsapiens, as.character(args[1]),start=as.numeric(args[6]), end=end(full.A)[[1]])
	c = getSeq(Hsapiens, as.character(args[1]),start=as.numeric(args[3]), end=as.numeric(args[4]))
	c.new = paste(rep(c,as.numeric(args[7])),collapse="")
	fileConn <- file(as.character(paste("Duplicated.",args[1],".genome.fasta",sep="")))
	sequence = as.character(paste(a,c.new,b,sep=""))
       	name = paste(">",args[1],sep="")
       	header = as.character(paste(name,sequence,sep="\n"))
      	writeLines(header,fileConn)
       	close(fileConn)
}
