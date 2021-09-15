args<-commandArgs(TRUE)
cor<-args[1]
input<-args[2]
out<-args[3]

library(fdrtool)
if(cor=='1'){
	data<-read.table(file=input,sep='\t',header=F)
	result<-fdrtool(data[,3], statistic="correlation", plot=FALSE, color.figure=FALSE)
	datanew<-cbind(data,result$qval,result$lfdr)
	write.table(datanew,file=out,sep='\t',quote=F,row.names=F,col.names=F)
}else{
	data<-read.table(file=input,sep='\t',header=F)
	result<-fdrtool(data[,7], statistic="pvalue", plot=FALSE, color.figure=FALSE)
	resultnew<-cbind(data,result$qval,result$lfdr)
	#resultout<-subset(resultnew,resultnew[,7]<1e-5&resultnew[,9]<0.001)
	write.table(resultnew,file=out,sep='\t',quote=F,row.names=F,col.names=F)
}




