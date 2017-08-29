library(gplots)
args <- commandArgs(TRUE)
mat <- as.matrix(read.table(args[1], header=FALSE))
feature <- read.table(args[2],header=FALSE)
x <- feature$V1
y <- feature$V1
z <- mat
z[z > quantile(z[z!=0],0.95)] <- min(z[z > quantile(z[z!=0],0.95)])
mx <- quantile(z[z!=0],1)
png(filename=args[3],width=8000,height=7000,res=1000)
filled.contour(x,y,z,zlim=range(1,mx), nlevels = 100, color.palette = colorRampPalette(c('red','gold2','yellow')), main="")
dev.off()
