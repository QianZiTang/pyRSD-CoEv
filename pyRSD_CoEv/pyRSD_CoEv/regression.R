args<-commandArgs(TRUE)
ratiofile<-args[1]
pcafile<-args[2]
output<-args[3]

data_evec<-read.table(file=pcafile,sep='\t',header=T,row.names = 1)
data_rate<-read.table(file=ratiofile,header=F,row.names = 1)

newdata_rate<-t(data_rate)
#newdata_rate[is.na(newdata_rate)]<-100
newdata_rate[is.na(newdata_rate)]<-max(data_rate,na.rm=TRUE)
new<-cbind(newdata_rate,data_evec)


getresi<-function(x){
  res<-resid(lm(x~PC1+PC2,data = new))
  #resall<-rbind(resall,res)
}

out<-apply(newdata_rate,2,getresi)  
out<-t(out)

write.table(out, file=output, sep='\t', quote=F, row.names=T, col.names=F)  
