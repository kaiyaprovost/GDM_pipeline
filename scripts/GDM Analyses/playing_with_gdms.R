df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate/bivariate_gdm_results.csv",
                header=T,sep="\t",fill=F)



correlation=cor(df[,c(
  #"number.fst.peaks","number.dxy.peaks",
  #        "Number.islands","number.sweeps","number.dxy.lows",
          "TAJIMAS","PERCENT_MISSING","CHROM_LENGTH",
          "MEAN_FST","MEAN_DXY","MEAN_RECOMB","prop.fst.peak",
          "prop.dxy.peak","prop.isl","prop.swp","prop.dxy.low")],
          use="pairwise.complete.obs")
corrplot::corrplot(correlation,method="ellipse",order="hclust",diag=F,
                   addrect=4)



y=(aggregate((df$number.sweeps+0.1) ~ df$MOSTA.1MODEL2,FUN=function(x){mean(x,na.rm=T)}))
z=aggregate((df$number.sweeps+0.1) ~ df$MOSTA.1MODEL2,FUN=function(x){sd(x,na.rm=T)})
barCenters=barplot(y[,2],ylim=c(0,100),names=y[,1],ylab="mean num sweeps")
segments(barCenters, y[,2]-z[,2], barCenters, y[,2]+z[,2], lwd=2)

y=(aggregate((df$number.sweeps+0.1) ~ df$MOSTA.1MODEL2,FUN=function(x){mean(log(x+0.1),na.rm=T)}))
barCenters=barplot(10^y[,2],names=y[,1],ylab="mean num sweeps (log calculated)")

y=(aggregate((df$prop.swp) ~ df$MOSTA.1MODEL2,FUN=function(x){mean(x,na.rm=T)}))
z=aggregate((df$prop.swp) ~ df$MOSTA.1MODEL2,FUN=function(x){sd(x,na.rm=T)})
barCenters=barplot(y[,2],names=y[,1],ylab="mean prop sweeps")


## PROPORTIONS
{
  ## not log 
  boxplot(df$prop.fst.peak~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
  vioplot::vioplot(df$prop.fst.peak~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))
  
  modelA=aov(df$prop.fst.peak[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
  summary(modelA) ## not sig
  TukeyHSD(modelA) 
  
  boxplot(df$prop.dxy.peak~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
  vioplot::vioplot(df$prop.dxy.peak~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))
  
  modelA=aov(df$prop.dxy.peak[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
  summary(modelA) ## sig
  TukeyHSD(modelA) ## mix-iba+ibe+ibh
  
  boxplot(df$prop.isl~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
  vioplot::vioplot(df$prop.isl~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))
  
  modelA=aov(df$prop.isl[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
  summary(modelA) ## not sig
  TukeyHSD(modelA) 
  
  boxplot(df$prop.swp~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
  vioplot::vioplot(df$prop.swp~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))
  
  modelA=aov(df$prop.swp[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
  summary(modelA) ## not sig
  TukeyHSD(modelA) 
  
  boxplot(df$prop.dxy.low~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
  vioplot::vioplot(df$prop.dxy.low~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))
  
  modelA=aov(df$prop.dxy.low[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
  summary(modelA) ## not sig
  TukeyHSD(modelA) 
}

## RAW NUMBERS 
{
## not log 
boxplot(df$number.fst.peaks~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(df$number.fst.peaks~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(df$number.fst.peaks[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(df$number.dxy.peaks~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(df$number.dxy.peaks~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(df$number.dxy.peaks[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(df$Number.islands~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(df$Number.islands~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(df$Number.islands[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(df$number.sweeps~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(df$number.sweeps~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(df$number.sweeps[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(df$number.dxy.lows~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(df$number.dxy.lows~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(df$number.dxy.lows[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

## log 
boxplot(log10(df$number.fst.peaks+0.1)~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.fst.peaks+0.1)~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.fst.peaks[df$MOSTA.1MODEL2!=""]+0.1) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$number.dxy.peaks+0.1)~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.dxy.peaks+0.1)~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.dxy.peaks[df$MOSTA.1MODEL2!=""]+0.1) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$Number.islands+0.1)~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$Number.islands+0.1)~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$Number.islands[df$MOSTA.1MODEL2!=""]+0.1) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$number.sweeps+0.1)~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.sweeps+0.1)~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.sweeps[df$MOSTA.1MODEL2!=""]+0.1) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$number.dxy.lows+0.1)~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.dxy.lows+0.1)~df$MOSTA.1MODEL2,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.dxy.lows[df$MOSTA.1MODEL2!=""]+0.1) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 
}

## chrom length
boxplot(log10(df$CHROM_LENGTH)~df$MOSTA.1MODEL2,las=2,varwidth=T,border="grey")
modelA=aov(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 
modelA=aov(log10(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 





## looking at differences across species
boxplot(log10(df$number.fst.peaks+0.1)~df$SPECIES,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.fst.peaks+0.1)~df$SPECIES,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.fst.peaks[df$MOSTA.1MODEL2!=""]+0.1) ~ df$SPECIES[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$number.dxy.peaks+0.1)~df$SPECIES,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.dxy.peaks+0.1)~df$SPECIES,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.dxy.peaks[df$MOSTA.1MODEL2!=""]+0.1) ~ df$SPECIES[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$Number.islands+0.1)~df$SPECIES,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$Number.islands+0.1)~df$SPECIES,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$Number.islands[df$MOSTA.1MODEL2!=""]+0.1) ~ df$SPECIES[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

boxplot(log10(df$number.sweeps+0.1)~df$SPECIES,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.sweeps+0.1)~df$SPECIES,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.sweeps[df$MOSTA.1MODEL2!=""]+0.1) ~ df$SPECIES[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

modelA1=aov(log10(df$number.sweeps[df$MOSTA.1MODEL2!=""]+0.1) ~ 
              df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""] + df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""])
summary(modelA1) ## not sig
TukeyHSD(modelA1) 

modelB=aov(log10(df$number.sweeps[df$MOSTA.1MODEL2!=""]+0.1) ~ 
             df$SPECIES[df$MOSTA.1MODEL2!=""] + df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelB) ## not sig
TukeyHSD(modelB) 

modelB1=aov(log10(df$number.sweeps[df$MOSTA.1MODEL2!=""]+0.1) ~ 
             df$SPECIES[df$MOSTA.1MODEL2!=""] + df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""] + df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""])
summary(modelB1) ## not sig
TukeyHSD(modelB1) 


boxplot(log10(df$number.dxy.lows+0.1)~df$SPECIES,las=2,varwidth=T,border="grey")
vioplot::vioplot(log10(df$number.dxy.lows+0.1)~df$SPECIES,las=2,colMed="red",add=T,col=rgb(0,0,0,0.1))

modelA=aov(log10(df$number.dxy.lows[df$MOSTA.1MODEL2!=""]+0.1) ~ df$SPECIES[df$MOSTA.1MODEL2!=""])
summary(modelA) ## not sig
TukeyHSD(modelA) 

