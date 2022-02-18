df=read.table("/Users/kprovost/gdm_weighted_mean_calc.txt",header=T,sep="\t")

for(spp in sort(unique(df$SPECIES))){
  for(col in 5:ncol(df)){
  temp = df[,c(1,2,3,col)]
  coln = colnames(df)[col]  
  
  wm=weighted.mean(df[df$SPECIES==spp,col],df$CHROM_LENGTH[df$SPECIES==spp],na.rm=T)
  m=mean(df[df$SPECIES==spp,col],na.rm=T)
  wstdv=sqrt(Hmisc::wtd.var(df[df$SPECIES==spp,col],df$CHROM_LENGTH[df$SPECIES==spp],na.rm=T))
  stdv=sd(df[df$SPECIES==spp,col],na.rm=T)
  
  print(paste(spp,coln,wm,"+",wstdv))
  }
}


