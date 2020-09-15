## testing ibd levels 

setwd("~")

genes=list.files(path="/Users/kprovost/Downloads/Dropbox/",
           pattern="_distancematrix_FULLGENOME.csv.converted",full.names = T)

geogs=list.files(path="/Users/kprovost/Downloads/Dropbox/",
           pattern="_distancematrix_geog.csv.converted",full.names = T)

for(i in 1:10){
  spp=strsplit(basename(genes[i]),"_")[[1]][1]
  print(spp)
  print("gene")
  gene = (read.table(genes[i],fill=F,header=T,row.names = 1,sep=","))
  print("geog")
  geog = (read.table(geogs[i],fill=F,header=T,row.names = 1,sep=","))
  
  
  genedf=NULL
  for (j in 1:ncol(gene)){
    toadd=cbind(rownames(gene),rep(colnames(gene)[j],nrow(gene)),as.numeric(gene[,j])) 
    genedf = rbind(genedf,toadd)
  }
  colnames(genedf) = c("row","col","gene")
  
  geogdf=NULL
  for (j in 1:ncol(geog)){
    toadd=cbind(rownames(geog),rep(colnames(geog)[j],nrow(geog)),as.numeric(geog[,j])) 
    geogdf = rbind(geogdf,toadd)
  }
  colnames(geogdf) = c("row","col","geog")
  
  png(paste(spp,"_gene_geog_ibd.png",sep=""))
  merged=(unique(merge(geogdf,genedf)))
  mod=lm(as.numeric(merged$gene)~as.numeric(merged$geog))
  pval=as.numeric(summary(mod)$coefficients[,4][2])
  dotted=2-as.numeric(pval<0.05)
  plot(as.numeric(merged$geog),as.numeric(merged$gene),ylab=signif(pval,3),
       xlab=signif(summary(mod)$adj.r.sq,2),main=spp)
  abline(mod,col="red",lty=dotted)
  dev.off()
  
}
