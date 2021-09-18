args=commandArgs(T)
mydata<-read.table(args[1],sep='\t',header=T,stringsAsFactors = F)
mydata2<-read.table(args[2],sep='\t',stringsAsFactors = F,header = F)


pdf(args[3],w=6,h=6)
d1<-density(mydata$TajimaD[mydata$TajimaD!="NaN"])
d2<-density(mydata2$V4[mydata2$V4!="NaN"])
x_min<-floor(min(range(d1$x)[1],range(d2$x)[1]))
x_max<-round(max(range(d1$x)[2],range(d2$x)[2]))

par(mar=c(5,5,4,5)+0.1)
plot(d2$x,d2$y,type = 'l',xlim = c(-1.2,3.5),col='red',xlab = "Taijima's D",ylab='Density (PASS regions)',las=1)
par(new=T)
plot(d1$x,d1$y,type = 'l',xlim = c(-1.2,3.5),ylim = c(0,2),col='black',xlab = "",ylab="",axes=F)
axis(4,las=1)
mtext("Density (Whole genome)",side=4,line=3)
legend("topright", inset=.05, c("PASS regions","Whole genome"),lty=c(1, 1), col=c("red", "black"),bty='n')
dev.off()
