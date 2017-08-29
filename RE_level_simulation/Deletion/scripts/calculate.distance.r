library(splitstackshape)
args <- commandArgs(trailingOnly = TRUE)
re.all = read.table(as.character(args[1]),header=F)
chr = read.table("../scripts/chr.list",header=F)
r = read.table("../scripts/IMR90.Distance.100Mb.prob",header=F)
imr90.dist = r$V1
imr90.prob = r$V2/sum(r$V2)
imr90.spline = smooth.spline(imr90.dist,imr90.prob,cv=T)
data = c()
i = as.numeric(args[3])
while (i <= as.numeric(args[3])){
	re = subset(re.all,V1==as.character(chr$V1[i]))
	pred.re = c()
	pred.count = c()
	master = c()
	master.chr = c()
	j = 1
	while (j <= length(re$V2)){
		pos = re$V2[j:length(re$V2)]
		d = abs(pos-pos[1])
		d[d > 100000000] = 100000000
		pred = 1-round(predict(imr90.spline,d)$y,5)
		rnb.count = rnbinom(n=length(pos),size=1500,prob=pred)
		pred.re = pos[which(rnb.count > 0)]
		pred.count = rnb.count[rnb.count > 0]
		master = rep(pos[1],length(pred.re))
		master.chr = rep(as.character(chr$V1[i]),length(pred.re))
		tmp = data.frame(master.chr,master,pred.re,pred.count)
		tmp = tmp[master!=pred.re,]
		tmp = expandRows(tmp, "pred.count")
		write.table(tmp,file=as.character(paste("IMR90.",chr$V1[i],".",args[2],".re.interaction.txt",sep="")),row.names=F,col.names=F,sep="\t",quote=F,append=TRUE)
		j = j+1
	}
	i = i+1
}
