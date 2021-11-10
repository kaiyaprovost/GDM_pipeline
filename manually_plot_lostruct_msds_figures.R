path = "/Users/kprovost/Downloads/Dropbox/"
files = list.files(path=path,pattern="coords.csv",full.names = T)
for(file in files){
  print(file)
  df = read.table(file,header=T,sep=",")
  spp = strsplit(basename(file),"\\.")[[1]][2]
  png(paste(spp,"_MSDS_plot_manual_lostruct.png",sep=""),height=3,width=3,units="in",res=600)
  par(mar=c(4,4,1,0))
  plot(df$MDS1,df$MDS2,col=df$ccols,main="",pch=as.numeric(as.factor(df$ccols)),ylab="MDS2",xlab="MDS1",las=2)
  dev.off()
  
  png(paste(spp,"_MSDS_boxplot_manual_lostruct.png",sep=""))
  par(mar=c(4,4,1,0),mfrow=(c(2,1)))
  boxplot(df$MDS1~df$chrom,main=spp,ylab="MDS1",xlab="Chrom",las=2)
  boxplot(df$MDS2~df$chrom,main=spp,ylab="MDS2",xlab="Chrom",las=2)
  dev.off()
  
  
}
