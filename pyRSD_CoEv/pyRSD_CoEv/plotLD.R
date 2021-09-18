args=commandArgs(T)
passfile<-read.table(args[1])
allgenome<-read.table(args[2])
outname<-args[3]
pdf(outname)
plot(data[,1]/1000,passfile[,2],type="l",col="red",main="LD decay",xlab="Distance(Kb)",xlim=c(0,500),ylim=c(0,1),ylab=expression(r^{2}),bty="n")
par(new=T)
plot(data2[,1]/1000,allgenome[,2],type="l",col="blue",main="",xlab="",xlim=c(0,500),ylim=c(0,1),ylab="",bty="n",axes=F)
dev.off()
