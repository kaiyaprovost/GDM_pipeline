library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
library(corrplot)
library(RColorBrewer)

morph <- read.table("~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_24february2018.csv",
                    sep=",",header=T,stringsAsFactors = F)

## notes on correcting for body mass:
## To correct for allometric scaling of beak size, we included
## body mass as a covariate in each analysis. This approach is
## preferred among contemporary phylogenetic comparative studies 
## because the use of residuals may cause collinearity issues
# (Freckleton 2009; Symonds and Tattersall 2010; Baab et al. 2014;
# Benson-Amram et al. 2015).

## Morphological variables and
# body mass were log transformed before analysis. Since
# nonnormality was not a problem with latitude, altitude,
# and Tmin, these were untransformed, except in the case of
# the group analyses of gulls and terns, where Nudds and
# Oswald’s log-transformed measure of latitude was used.
# To control for body mass, we calculated multiple regres-
# sions in COMPARE, predicting bill length (or leg-element
# length), with both latitude (or altitude or minimum tem-
# perature) and body mass as predictors.

#nobad = morph[morph$CONDITION=="",]
#noblank1 = nobad[nobad$BEAKAREAVOLUME!="",]
#noblank2 = noblank1[!(is.na(noblank1$BEAKAREAVOLUME)),]
#noblank = noblank2[!(is.na(noblank2$KIPPSINDEX)),]
#males = noblank[noblank$SEX=="MALE",]

agg = aggregate(cbind(BILL.HEIGHT,BILL.LENGTH,BILL.WIDTH,TARSUS.LENGTH,WING.LENGTH.TO.PRIMARIES,
                      WING.LENGTH.TO.SECONDARIES,BAD.MEASUREMENT,LAT,LONG,KIPPSINDEX,STRUCTURE,BEAKBASEAREA,BEAKVOL,BEAKLATERALSURFACE,
                      BEAKTOTALSURFACE,BEAKAREAVOLUME) ~ CATALOG.NUMBER, data = morph, FUN=mean,na.action=na.pass)
agg2 = aggregate(cbind(SPECIES,WHICH.SIDE.OF.CFB,COUNTRY,STATE,COUNTY,LOCALITY,SEX,AGE,CONDITION,MEASURER,DATE,NOTES,
                       GENUS,SPP,SUBSPP,GEOREF.BY,LOCALITY2)~CATALOG.NUMBER,data=morph,
                 FUN=function(x) {paste(unique(x),collapse="; ")},na.action=na.pass)
agg3 = aggregate(cbind(MEASUREMENT)~CATALOG.NUMBER,data=morph,FUN=max,na.action=na.pass)

mergeAgg = merge(agg3,agg,by="CATALOG.NUMBER",all=T)
mergeAgg = merge(mergeAgg,agg2,by="CATALOG.NUMBER",all=T)

kipps = function(primary,secondary) {
  #kipps index calculated (primary_length-secondary_length/primary_length*100)
  return(((primary-secondary)/primary)*100)
}
beakbase = function(height,width){
  #beak base area calculated as area base cone (pi*beak_height*beak_width)
  return(pi*height*width)
}
beakvolume = function(height,width,length){
  #beak volume calculated as beak vol cone (1/3*beak_base*beak_length)
  base = beakbase(height,width)
  return((1/3)*base*length)
}
beakappxlatsurf = function(height,width,length){
  #beak lateral surface area calculated as appx lateral surf area cone ((beak_width+beak_height)/4)*beak_length*pi
  return(((width+height)/4)*length*pi)
}
beaktotsurf = function(height,width,length){
  #beak total surface is beak_lateral_surface + beak_base_area
  base = beakbase(height,width)
  lat = beakappxlatsurf(height,width,length)
  return(base+lat)
}
beaklatSAVR = function(height,width,length){
  #beak area volume is surface area / beak volume
  vol = beakvolume(height,width,length)
  lat = beakappxlatsurf(height,width,length)
  return(lat/vol)
}
beaktotSAVR = function(height,width,length){
  #beak area volume is surface area / beak volume
  vol = beakvolume(height,width,length)
  tot = beaktotsurf(height,width,length)
  return(tot/vol)
}
beakgeom = function(height,width,length){
  base = beakbase(height,width)
  vol = beakvolume(height,width,length)
  lat = beakappxlat(height,width,length)
  tot = beaktotsurf(height,width,length)
  latvolratio = beaklatSAVR(height,width,length)
  totvolratio = beaktotSAVR(height,width,length)
  x = as.data.frame(cbind(base,vol,lat,tot,latvolratio,totvolratio))
  return(x)
}

kipp = kipps(primary=mergeAgg$WING.LENGTH.TO.PRIMARIES,secondary=mergeAgg$WING.LENGTH.TO.SECONDARIES)
beaks = beakgeom(height=mergeAgg$BILL.HEIGHT,width=mergeAgg$BILL.WIDTH,length=mergeAgg$BILL.LENGTH)
calculated = cbind(kipp,beaks,mergeAgg$CATALOG.NUMBER)
names(calculated)
names(calculated) = c("KIPPSINDEX_CALC","BEAKBASEAREA_CALC","BEAKVOL_CALC",
                      "BEAKLATERALSURFACE_CALC","BEAKTOTALSURFACE_CALC",
                      "BEAKAREAVOLUME_CALC","BEAKTOTAREAVOLUME_CALC","CATALOG.NUMBER")
names(calculated)
mergeAgg = merge(mergeAgg,calculated,by="CATALOG.NUMBER",all=T)

write.csv(mergeAgg,file=paste("~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_",format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
          na="")

morph = read.table("~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_16_March_2018.csv",
                   sep=",",header=T)

####

## control for body size

log = log(morph[,c(4:9,37:43)])
names(log) = sapply(names(log),simplify=T,USE.NAMES = F, FUN=function(x) {paste("LOG",x,sep="_")})
log = cbind(morph[,c(1:3,10:12,14,20:36)],log)

write.csv(log,file=paste("~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_LOG_",format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
          na="")

model = lm(morph$BILL.HEIGHT~morph$TARSUS.LENGTH)
plot(morph$TARSUS.LENGTH,morph$BILL.HEIGHT)
abline(model)

model = lm(log$LOG_BILL.HEIGHT~log$LOG_TARSUS.LENGTH)
plot(log$LOG_TARSUS.LENGTH,log$LOG_BILL.HEIGHT)
abline(model)

## get residuals after log

model$residuals

res_table = log[,2:24]
rownames(res_table) = log$X
for (i in 25:length(names(log))) {
  name = names(log)[i]
  print(name)
  model = lm(log[,i]~log$LOG_TARSUS.LENGTH)
  res = model$residuals
  res = unname(res)
  x = merge(res_table,res,by="row.names",all.x=TRUE)
  names(x)[i] = paste("RES_",name,sep="")
  rownames(x) = x$Row.names
  res_table = x[2:length(colnames(x))]
}

write.csv(res_table,
          file=paste("~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_LOG_residuals_",
                     format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
          na="")

## check to see if residuals are correlated
corrplot(cor(res_table[,24:36],use="pairwise.complete.obs"),
         method="ellipse",type="upper",diag=T)



## calculate residuals before log and compare beak/wing
mor1 = morph[,c(4:9,37:43)]
mor1 = cbind(morph[,c(1:3,10:12,14,20:36)],mor1)


nolog_res = mor1[,2:24]
rownames(nolog_res) = mor1$X
for (i in 25:length(names(mor1))) {
  name = names(mor1)[i]
  print(name)
  model = lm(mor1[,i]~mor1$TARSUS.LENGTH)
  res = model$residuals
  res = unname(res)
  x = merge(nolog_res,res,by="row.names",all.x=TRUE)
  names(x)[i] = paste("RES_",name,sep="")
  rownames(x) = x$Row.names
  nolog_res = x[2:length(colnames(x))]
}

write.csv(nolog_res,
          file=paste("~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphomor1y_AGGREGATED_residuals_",
                     format(Sys.time(),"%d_%B_%Y"),".csv",sep=""),
          na="")

corrplot(cor(nolog_res[,24:36],use="pairwise.complete.obs"),
         method="ellipse",type="upper",diag=T)

nobad = nolog_res[nolog_res$CONDITION=="",]
#noblank1 = nobad[nobad$RES_BEAKAREAVOLUME_CALC!="",]
#noblank2 = noblank1[!(is.na(noblank1$RES_BEAKAREAVOLUME_CALC)),]
#noblank = noblank2[!(is.na(noblank2$RES_KIPPSINDEX_CALC)),]
#males = noblank[noblank$SEX=="MALE",]
males = nobad[nobad$SEX=="MALE",]

#males$SPP<-factor(males$SPP, levels=c("CARDINALIS","SINUATUS","FUSCA",
#                                      "CURVIROSTRE","CRISSALE/DORSALE",
#                                      "BILINEATA","BRUNNEICAPILLUS","NITENS",
#                                      "FLAVICEPS","BELLII","MELANURA"))

males$SPP = factor(males$SPP, levels=c("SINUATUS","CARDINALIS","FUSCA",
                                       "BILINEATA","NITENS","FLAVICEPS",
                                       "MELANURA","CRISSALE/DORSALE",
                                       "CURVIROSTRE","BRUNNEICAPILLUS","BELLII"))

names=c("SIN","CAR","FUS","BIL","NIT","FLA","MEL","CRI","CUR","BRU","BEL")
col=c("white","cyan","cyan","white","white","red","grey","white","cyan","red","cyan")


png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume_RESIDUALS.png")
boxplot(as.numeric(males$RES_BEAKAREAVOLUME_CALC)~males$SPP,horizontal=F,
        names=names,las=2,
        xlab="Spp",ylab="Corr Lat Surf Area (mm2) / Corr Beak Vol (mm3)",
        col=col)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfArea_RESIDUALS.png")
boxplot(as.numeric(males$RES_BEAKLATERALSURFACE_CALC)~males$SPP,horizontal=F,
        names=names,las=2,xlab="Spp",ylab="Corrected Lat Surface Area (mm2)",
        col=col)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume_RESIDUALS.png")
boxplot(as.numeric(males$RES_BEAKVOL_CALC)~males$SPP,horizontal=F,
        names=names,las=2,xlab="Spp",ylab="Corrected Log Beak Vol (mm3)",
        col=col)
dev.off()

summary(males$RES_BEAKAREAVOLUME_CALC)
meanSAV = aggregate(males$RES_BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=mean)
medSAV =  aggregate(males$RES_BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=median)

chi = males[males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]
son = males[males$WHICH.SIDE.OF.CFB=="SONORAN",]
unk = males[males$WHICH.SIDE.OF.CFB=="UNCLEAR",]
no_unk = males[males$WHICH.SIDE.OF.CFB!="UNCLEAR",]

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfArea_overlap_byDesert_RESIDUALS.png")
boxplot(son$RES_BEAKLATERALSURFACE_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-100,200),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_BEAKLATERALSURFACE_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-100,200),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume_overlap_byDesert_RESIDUALS.png")
boxplot(son$RES_BEAKVOL_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-800,2100),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_BEAKVOL_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-800,2100),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume_overlap_byDesert_RESIDUALS.png")
boxplot(son$RES_BEAKAREAVOLUME_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-0.25,0.5),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_BEAKAREAVOLUME_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-0.25,0.5),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_kipps_overlap_byDesert_RESIDUALS.png")
boxplot(son$RES_KIPPSINDEX_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-7,11),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_KIPPSINDEX_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-7,11),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

for (spp in unique(males$SPP)) {
  for (des in unique(males$WHICH.SIDE.OF.CFB)) {
    #print(paste(spp,des))
    tab = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB==des,c(24:36)]
    m = as.data.frame(sapply(tab,mean))
    #s = sapply(tab,sd)
    #names(spp) = "SPP"
    #names(des) = "WHICH.SIDE.OF.CFB"
    sppdes = paste(substring(spp,1,3),substring(des,1,3),sep="_")
    names(m) = sppdes
    tm = t(m)
    write(paste("SPP_DES",names(tab),collapse=",",sep=","),"~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_RESIDUALS_MEAN.csv",
          append=T,sep=",")
    write.table(tm,"~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_RESIDUALS_MEAN.csv",
                append=T,sep=",",col.names=F)
  }
}


for (spp in unique(males$SPP)) {
  print(spp)
  write(spp,"/Users/kprovost/Documents/Classes/Stephanie/output_RES.txt",append=T,sep="\t",ncolumns=3)
  tabS = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="SONORAN",c(24:36)]
  tabC = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",c(24:36)]
  for (colnum in 1:length(colnames(tabS))) {
    write(colnames(tabS)[colnum],"/Users/kprovost/Documents/Classes/Stephanie/output_RES.txt",append=T,sep="\t",ncolumns=3)
    print(colnum)
    if ((nrow(tabC) != 0) && (nrow(tabS) != 0) && sum(!is.na(tabC[,colnum])) != 0 && sum(!is.na(tabS[,colnum]))) {
      test = t.test(tabS[,colnum],tabC[,colnum])
      write(test$p.value,"/Users/kprovost/Documents/Classes/Stephanie/output_RES.txt",append=T,sep="\n",ncolumns=3)
    } else {
      write("#####","/Users/kprovost/Documents/Classes/Stephanie/output_RES.txt",append=T,sep="\n",ncolumns=3)
    }
  }
}

testtable = read.table("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test_RES.txt",header=T,
                       sep="\t")



##### 

## permute son vs chi values and test 
## start off with melanura and beak area volume
mel1 = mor1[mor1$SPP=="MELANURA",]
bav = (mel1$BEAKAREAVOLUME_CALC)
des = mel1$WHICH.SIDE.OF.CFB

N=5
x = replicate(N,sample(des),simplify=T)
y = as.data.frame(cbind(as.numeric(bav),as.character(des),x),stringsAsFactors = F)
colnames(y) = c("VALUE","ORIGINAL",paste("REP",1:N,sep="_"))
y$VALUE = as.numeric(y$VALUE)
t.test(y$VALUE[y$ORIGINAL=="SONORAN"],y$VALUE[y$ORIGINAL=="CHIHUAHUAN"])


#####






## comparing beak and wing stuff after fix for tarsus (LOG)

nobad = res_table[res_table$CONDITION=="",]
#noblank1 = nobad[nobad$RES_LOG_BEAKAREAVOLUME_CALC!="",]
#noblank2 = noblank1[!(is.na(noblank1$RES_LOG_BEAKAREAVOLUME_CALC)),]
#noblank = noblank2[!(is.na(noblank2$RES_LOG_KIPPSINDEX_CALC)),]
#males = noblank[noblank$SEX=="MALE",]
males = nobad[nobad$SEX=="MALE",]

#males$SPP<-factor(males$SPP, levels=c("CARDINALIS","SINUATUS","FUSCA",
#                                      "CURVIROSTRE","CRISSALE/DORSALE",
#                                      "BILINEATA","BRUNNEICAPILLUS","NITENS",
#                                      "FLAVICEPS","BELLII","MELANURA"))

males$SPP = factor(males$SPP, levels=c("SINUATUS","CARDINALIS","FUSCA",
                                       "BILINEATA","NITENS","FLAVICEPS",
                                       "MELANURA","CRISSALE/DORSALE",
                                       "CURVIROSTRE","BRUNNEICAPILLUS","BELLII"))

names=c("SIN","CAR","FUS","BIL","NIT","FLA","MEL","CRI","CUR","BRU","BEL")
col=c("white","cyan","cyan","white","white","red","grey","white","cyan","red","cyan")


png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume_LOG_RESIDUALS.png")
boxplot(as.numeric(males$RES_LOG_BEAKAREAVOLUME_CALC)~males$SPP,horizontal=F,
        names=names,las=2,
        xlab="Spp",ylab="Corr Log Lat Surf Area (mm2) / Corr Log Beak Vol (mm3)",
        col=col)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfArea_LOG_RESIDUALS.png")
boxplot(as.numeric(males$RES_LOG_BEAKLATERALSURFACE_CALC)~males$SPP,horizontal=F,
        names=names,las=2,xlab="Spp",ylab="Corrected Log Lat Surface Area (mm2)",
        col=col)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume_LOG_RESIDUALS.png")
boxplot(as.numeric(males$RES_LOG_BEAKVOL_CALC)~males$SPP,horizontal=F,
        names=names,las=2,xlab="Spp",ylab="Corrected Log Beak Vol (mm3)",
        col=col)
dev.off()

summary(males$RES_LOG_BEAKAREAVOLUME_CALC)
meanSAV = aggregate(males$RES_LOG_BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=mean)
medSAV =  aggregate(males$RES_LOG_BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=median)

chi = males[males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]
son = males[males$WHICH.SIDE.OF.CFB=="SONORAN",]
unk = males[males$WHICH.SIDE.OF.CFB=="UNCLEAR",]
no_unk = males[males$WHICH.SIDE.OF.CFB!="UNCLEAR",]

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfArea_overlap_byDesert_LOG_RESIDUALS.png")
boxplot(son$RES_LOG_BEAKLATERALSURFACE_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-1,1.25),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_LOG_BEAKLATERALSURFACE_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-1,1.25),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume_overlap_byDesert_LOG_RESIDUALS.png")
boxplot(son$RES_LOG_BEAKVOL_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-2,2),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_LOG_BEAKVOL_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-2,2),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume_overlap_byDesert_LOG_RESIDUALS.png")
boxplot(son$RES_LOG_BEAKAREAVOLUME_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-1,1),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_LOG_BEAKAREAVOLUME_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-1,1),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_kipps_overlap_byDesert_LOG_RESIDUALS.png")
boxplot(son$RES_LOG_KIPPSINDEX_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(-1,1),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1,names=names)
boxplot(chi$RES_LOG_KIPPSINDEX_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(-1,1),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5,names=names)
dev.off()

for (spp in unique(males$SPP)) {
  for (des in unique(males$WHICH.SIDE.OF.CFB)) {
    #print(paste(spp,des))
    tab = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB==des,c(24:36)]
    m = as.data.frame(sapply(tab,mean))
    #s = sapply(tab,sd)
    #names(spp) = "SPP"
    #names(des) = "WHICH.SIDE.OF.CFB"
    sppdes = paste(substring(spp,1,3),substring(des,1,3),sep="_")
    names(m) = sppdes
    tm = t(m)
    write(paste("SPP_DES",names(tab),collapse=",",sep=","),"~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_LOG_RESIDUALS_MEAN.csv",
                append=T,sep=",")
    write.table(tm,"~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_LOG_RESIDUALS_MEAN.csv",
                append=T,sep=",",col.names=F)
  }
}


for (spp in unique(males$SPP)) {
  print(spp)
  write(spp,"/Users/kprovost/Documents/Classes/Stephanie/output_LOG_RES.txt",append=T,sep="\t",ncolumns=3)
  tabS = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="SONORAN",c(24:36)]
  tabC = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",c(24:36)]
  for (colnum in 1:length(colnames(tabS))) {
    write(colnames(tabS)[colnum],"/Users/kprovost/Documents/Classes/Stephanie/output_LOG_RES.txt",append=T,sep="\t",ncolumns=3)
    print(colnum)
    if ((nrow(tabC) != 0) && (nrow(tabS) != 0) && sum(!is.na(tabC[,colnum])) != 0 && sum(!is.na(tabS[,colnum]))) {
      test = t.test(tabS[,colnum],tabC[,colnum])
      write(test$p.value,"/Users/kprovost/Documents/Classes/Stephanie/output_LOG_RES.txt",append=T,sep="\n",ncolumns=3)
    } else {
      write("#####","/Users/kprovost/Documents/Classes/Stephanie/output_LOG_RES.txt",append=T,sep="\n",ncolumns=3)
    }
  }
}

testtable = read.table("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test_LOG_RES.txt",header=T,
                       sep="\t")







## comparing beak area volume etc -- check if working -- not fixed for tarsus
#####

nobad = morph[morph$CONDITION=="",]
noblank1 = nobad[nobad$BEAKAREAVOLUME_CALC!="",]
noblank2 = noblank1[!(is.na(noblank1$BEAKAREAVOLUME_CALC)),]
noblank = noblank2[!(is.na(noblank2$KIPPSINDEX_CALC)),]
males = noblank[noblank$SEX=="MALE",]

males$SPP<-factor(males$SPP, levels=c("CARDINALIS","SINUATUS","FUSCA","CURVIROSTRE","CRISSALE/DORSALE","BILINEATA","BRUNNEICAPILLUS","NITENS","FLAVICEPS","BELLII","MELANURA"))

png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume.png")
boxplot(as.numeric(males$BEAKAREAVOLUME_CALC)~males$SPP,horizontal=F,
        names=c("car","sin","fus","cur","cri","bil","bru","nit","fla","bel","mel"),
        xlab="Spp",ylab="Lateral Surface Area (mm2) / Beak Volume (mm3)",
        col=c("cyan","white","cyan","cyan","white","white","red","white","red","cyan","grey"))
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfArea.png")
boxplot(as.numeric(males$BEAKLATERALSURFACE_CALC)~males$SPP,horizontal=F,
        names=c("car","sin","fus","cur","cri","bil","bru","nit","fla","bel","mel"),
        col=c("cyan","white","cyan","cyan","white","white","red","white","red","cyan","grey"))
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume.png")
boxplot(as.numeric(males$BEAKVOL_CALC)~males$SPP,horizontal=F,
        names=c("car","sin","fus","cur","cri","bil","bru","nit","fla","bel","mel"),
        col=c("cyan","white","cyan","cyan","white","white","red","white","red","cyan","grey"))
dev.off()

summary(males$BEAKAREAVOLUME_CALC)
meanSAV = aggregate(males$BEAKAREAVOLUME_CALC, by=list(males$SPP), FUN=mean)

#####
#comparing between deserts

chi = males[males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",]
son = males[males$WHICH.SIDE.OF.CFB=="SONORAN",]
unk = males[males$WHICH.SIDE.OF.CFB=="UNCLEAR",]

no_unk = males[males$WHICH.SIDE.OF.CFB!="UNCLEAR",]

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfaceArea_byDesert.png")
boxplot(no_unk$BEAKLATERALSURFACE_CALC~no_unk$WHICH.SIDE.OF.CFB+no_unk$SPP,
        las=2,col=c("blue","green","grey"),
        names=c("car.C","car.S","","sin.C","sin.S","",
                "fus.C","fus.S","","cur.C","cur.S","",
                "cri.C","cri.S","","bil.C","bil.S","",
                "bru.C","bru.S","","nit.C","nit.S","",
                "fla.C","fla.S","","bel.C","bel.S","",
                "mel.C","mel.S",""),
        xlab="Species-desert",ylab="Beak lateral surface area")
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume_byDesert.png")
boxplot(no_unk$BEAKAREAVOLUME_CALC~no_unk$WHICH.SIDE.OF.CFB+no_unk$SPP,
        las=2,col=c("blue","green","grey"),
        names=c("car.C","car.S","","sin.C","sin.S","",
                "fus.C","fus.S","","cur.C","cur.S","",
                "cri.C","cri.S","","bil.C","bil.S","",
                "bru.C","bru.S","","nit.C","nit.S","",
                "fla.C","fla.S","","bel.C","bel.S","",
                "mel.C","mel.S",""),
        xlab="Species-desert",ylab="Beak surface area / volume")
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume_overlap_byDesert.png")
boxplot(no_unk$BEAKVOL_CALC~no_unk$WHICH.SIDE.OF.CFB+no_unk$SPP,
        las=2,col=c("blue","green","grey"),
        names=c("car.C","car.S","","sin.C","sin.S","",
                "fus.C","fus.S","","cur.C","cur.S","",
                "cri.C","cri.S","","bil.C","bil.S","",
                "bru.C","bru.S","","nit.C","nit.S","",
                "fla.C","fla.S","","bel.C","bel.S","",
                "mel.C","mel.S",""),
        xlab="Species-desert",ylab="Beak volume")
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakLateralSurfArea_overlap_byDesert.png")
boxplot(son$BEAKLATERALSURFACE_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(0,450),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1)
boxplot(chi$BEAKLATERALSURFACE_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(0,450),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakVolume_overlap_byDesert.png")
boxplot(son$BEAKVOL_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(0,3500),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1)
boxplot(chi$BEAKVOL_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(0,3500),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5)
dev.off()

png("~/Documents/Classes/Stephanie/aggregated_beakAreaVolume_overlap_byDesert.png")
boxplot(son$BEAKAREAVOLUME_CALC~son$SPP,border=rgb(0,1,0,0.60),ylim=c(0,0.7),add=F,col=rgb(0,1,0,0.30),
        las=2,boxwex=1)
boxplot(chi$BEAKAREAVOLUME_CALC~chi$SPP,border=rgb(0,0,1,0.60),ylim=c(0,0.7),add=T,col=rgb(0,0,1,0.30),
        las=2,boxwex=0.5)
dev.off()

for (spp in unique(males$SPP)) {
  for (des in unique(males$WHICH.SIDE.OF.CFB)) {
    #print(paste(spp,des))
    tab = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB==des,c(4:9,13,15:19,37:43)]
    m = as.data.frame(sapply(tab,mean))
    #s = sapply(tab,sd)
    #names(spp) = "SPP"
    #names(des) = "WHICH.SIDE.OF.CFB"
    sppdes = paste(substring(spp,1,3),substring(des,1,3),sep="_")
    names(m) = sppdes
    tm = t(m)
    write.table(tm,"~/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_TEST.csv",
                append=T,sep=",",col.names=F)
    }
  }

means_values = read.table("/Users/kprovost/Documents/Classes/Stephanie/KLP_Master_Spreadsheet_Morphology_AGGREGATED_means.csv",
                          sep="\t",header=T)
names(means_values)[1] = "SPPDES"
#par(ask=T)
for (i in 2:length(names(means_values))){
  filename = paste("/Users/kprovost/Documents/Classes/Stephanie/barplot_",
                   colnames(means_values[i]),".png",sep="")
  png(filename)
  barplot(means_values[,i],names=means_values[,1],las=2,
          col=c("green","blue","green","blue",
                "green","blue","green","blue",
                "green","blue","green","blue",
                "grey","green","blue",
                "green","blue","green","blue",
                "green","grey","green","blue"),
          space=c(0,0,0.3,0,0.3,0,0.3,0,
                  0.3,0,0.3,0,0,0.3,0,
                  0.3,0,0.3,0,0.3,0,0.3,0),
          ylab=colnames(means_values)[i])
  dev.off()
}
#par(ask=F)

for (spp in unique(males$SPP)) {
    write(spp,"/Users/kprovost/Documents/Classes/Stephanie/output.txt",append=T,sep="\t",ncolumns=3)
    tabS = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="SONORAN",c(4:9,13,15:19,37:43)]
    tabC = males[males$SPP==spp & males$WHICH.SIDE.OF.CFB=="CHIHUAHUAN",c(4:9,13,15:19,37:43)]
    for (colnum in 1:length(colnames(tabS))) {
      write(colnames(tabS)[colnum],"/Users/kprovost/Documents/Classes/Stephanie/output.txt",append=T,sep="\t",ncolumns=3)
      if ((nrow(tabC) != 0) && (nrow(tabS) != 0)) {
        test = t.test(tabS[,colnum],tabC[,colnum])
        write(test$p.value,"/Users/kprovost/Documents/Classes/Stephanie/output.txt",append=T,sep="\n",ncolumns=3)
      } else {
        write("#####","/Users/kprovost/Documents/Classes/Stephanie/output.txt",append=T,sep="\n",ncolumns=3)
      }
    }
}

testtable = read.table("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test.txt",header=T,
                       sep="\t")

png("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test.png")
ramp = brewer.pal(12,"Paired")
ramp[11] = "grey"
rain = rainbow(14)
cols = c("black","black",ramp)
pchs=c(0,21,22,23,24,25,0,1,2,5,6,11,3,4)
par(xpd = T, mar = par()$mar + c(0,0,0,3))
for (colnum in 2:length(names(testtable))){
  if (colnum == 2){
    plot(1:10,testtable[,colnum],xaxt="n",xlab="",
         ylab="p-value Son vs Chi t-test",
         ylim=c(0,1),pch=colnum-2,bg=cols[colnum],lwd=2)
    axis(1,at=1:10,labels=testtable$SPECIES,las=2)
    } else {
      points((1:10+((colnum-2)/25)-0.25),testtable[,colnum],col=cols[colnum],
             ylim=c(0,1),pch=colnum-2,bg=cols[colnum],lwd=2)
      }
  #boxplot(testtable[,colnum]~testtable$SPECIES,
  #        las=2,horizontal=T,add=add,border=cols[colnum],
  #        col=cols[colnum])
}
abline(h=0.05,col="red",lty=2,xpd=F)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),
       col="gray90",xpd=F)
legend(x=10.5,y=0.8,legend=c("BLSA:V","BBA","BLSA",
                        "BTSA","BTSA:V","BV","BH",
                        "BL","BW","KI","TL","WLP","WLS"),
       xpd=T,col=cols[2:14],bty="n",
       pch=0:12,pt.lwd=2)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()



png("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test_panel.png")
ramp = brewer.pal(12,"Paired")
ramp[11] = "grey"
rain = rainbow(14)
cols = c("black","black",ramp)
pchs=c(0,21,22,23,24,25,0,1,2,5,6,11,3,4)
names(testtable) = c("SPECIES","BLSA:V","BBA","BLSA",
                     "BTSA","BTSA:V","BV","BH",
                     "BL","BW","KI","TL","WLP","WLS")
par(mfrow=c(4,4),mar=c(3,4,0,0))
for (colnum in 2:length(names(testtable))){
    plot(1:10,testtable[,colnum],xaxt="n",xlab="",
         ylab=names(testtable[colnum]),
         ylim=c(0,1),pch=colnum-2,bg=cols[colnum],lwd=2)
    axis(1,at=1:10,labels=substr((testtable$SPECIES),1,3),las=2)
    # "red","blue","blue","blue","red",
    #"white","blue","white","white","white"
    abline(h=0.05,col="red",lty=2,xpd=F)
}
dev.off()


png("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test_zoom.png")
ramp = brewer.pal(12,"Paired")
ramp[11] = "grey"
rain = rainbow(14)
cols = c("black","black",ramp)
pchs=c(0,21,22,23,24,25,0,1,2,5,6,11,3,4)
par(xpd = T, mar = par()$mar + c(0,0,0,3))
for (colnum in 2:length(names(testtable))){
  if (colnum == 2){
    plot(1:10,testtable[,colnum],xaxt="n",xlab="",
         ylab="p-value Son vs Chi t-test",
         ylim=c(0,0.06),pch=colnum-2,bg=cols[colnum],lwd=2,xpd=F)
    axis(1,at=1:10,labels=testtable$SPECIES,las=2)
  } else {
    points((1:10+((colnum-2)/25)-0.25),testtable[,colnum],col=cols[colnum],
           ylim=c(0,0.06),pch=colnum-2,bg=cols[colnum],lwd=2,xpd=F)
  }
  #boxplot(testtable[,colnum]~testtable$SPECIES,
  #        las=2,horizontal=T,add=add,border=cols[colnum],
  #        col=cols[colnum])
}
abline(h=0.05,col="red",lty=2,xpd=F)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),
       col="gray90",xpd=F)
legend(x=10.5,y=0.05,legend=c("BLSA:V","BBA","BLSA",
                             "BTSA","BTSA:V","BV","BH",
                             "BL","BW","KI","TL","WLP","WLS"),
       xpd=T,col=cols[2:14],bty="n",
       pch=0:12,pt.lwd=2)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


png("/Users/kprovost/Documents/Classes/Stephanie/sonVsChit-test_zoom01.png")
ramp = brewer.pal(12,"Paired")
ramp[11] = "grey"
rain = rainbow(14)
cols = c("black","black",ramp)
pchs=c(0,21,22,23,24,25,0,1,2,5,6,11,3,4)
par(xpd = T, mar = par()$mar + c(0,0,0,3))
for (colnum in 2:length(names(testtable))){
  if (colnum == 2){
    plot(1:10,testtable[,colnum],xaxt="n",xlab="",
         ylab="p-value Son vs Chi t-test",
         ylim=c(0,0.01),pch=colnum-2,bg=cols[colnum],lwd=2,xpd=F)
    axis(1,at=1:10,labels=testtable$SPECIES,las=2)
  } else {
    points((1:10+((colnum-2)/25)-0.25),testtable[,colnum],col=cols[colnum],
           ylim=c(0,0.01),pch=colnum-2,bg=cols[colnum],lwd=2,xpd=F)
  }
  #boxplot(testtable[,colnum]~testtable$SPECIES,
  #        las=2,horizontal=T,add=add,border=cols[colnum],
  #        col=cols[colnum])
}
abline(h=0.05,col="red",lty=2,xpd=F)
abline(v=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),
       col="gray90",xpd=F)
legend(x=10.5,y=0.008,legend=c("BLSA:V","BBA","BLSA",
                              "BTSA","BTSA:V","BV","BH",
                              "BL","BW","KI","TL","WLP","WLS"),
       xpd=T,col=cols[2:14],bty="n",
       pch=0:12,pt.lwd=2)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()




#####

bellii = morph[morph$SPP=="bellii",] # done
bilineata = morph[morph$SPP=="bilineata",]
brunneicapillus = morph[morph$SPP=="brunneicapillus",] # done
cardinalis = morph[morph$SPP=="cardinalis",] #d one
chlorurus = morph[morph$SPP=="chlorurus",] #d one
crissalis = morph[morph$SPP=="crissalis",] #d one
curvirostre = morph[morph$SPP=="curvirostre",] # done
dorsale = morph[morph$SPP=="dorsale",] #d one
flaviceps = morph[morph$SPP=="flaviceps",] # done
fusca = morph[morph$SPP=="fusca",] # done
melanura = morph[morph$SPP=="melanura",] # done
nitens = morph[morph$SPP=="nitens",]
sinuatus = morph[morph$SPP=="sinuatus",] # done


boxplot(morph$kippsindex~morph$SPP,horizontal=F,
        names=c("0-bel","?-bil","1-bru","0-car","","","0-cur",
                "1-fla","0-fus","0-mel","?-nit","?-sin"))


morph$STATE[morph$STATE=="ARIZONA"] = "Arizona"
morph$STATE[morph$STATE=="TEXAS"] = "Texas"
morph$AGE[morph$AGE=="Juv."] = "juv"
morph = morph[!is.na(morph$SPP),]
morph = morph[!(morph$SPP=="chlorurus"),]
morph = morph[!(morph$SPP=="crissalis"),]
morph = morph[!(morph$SPP=="plumbea"),]
morph = morph[!(morph$SEX=="Female"),]
morph = morph[!(morph$AGE=="juv"),]
morph = morph[!is.na(morph$BILL.LENGTH..mm.),]
morph = morph[!is.na(morph$WING.LENGTH.TO.PRIMARIES..mm.),]
morph = morph[!is.na(morph$BILL.WIDTH..mm.),]
morph = morph[!is.na(morph$WING.LENGTH.TO.SECONDARIES..mm.),]
morph = morph[!is.na(morph$TARSUS.LENGTH..mm.),]
morph_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                     TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.,
                   data=morph,center=T,scale=T,na.action = "na.omit")

summary(morph_pca)
morph_with_pca = merge(morph,morph_pca$x,by="row.names")
corrplot(cor(morph_with_pca[c(11:16,29:33)]),
         method="number",type="upper",diag=F)

plot(morph_with_pca$PC1,morph_with_pca$PC3)
plot(morph_with_pca$PC2,morph_with_pca$PC3)


g <- ggbiplot(morph_pca, obs.scale = 1, var.scale = 1, 
              groups = morph$SPP, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
{
g = g +
  #scale_color_manual(name="", values=c("green", "yellow","pink","lightgreen","cyan",
  #                                     "red","darkgreen","magenta","orange","white","purple")) +
  #scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
  #                                     "red","cyan","red","white","white","white")) +
  scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
                                       "cyan","red","cyan","cyan","red","white")) +
   scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(aes(colour=morph$SPP, shape=morph$SPP))
g = g + ggtitle("")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/morph_pca_spp_enm_not_26feb2017.png",
    width=12,height=6,pointsize=1,units="cm",
    res=500)
#par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()



## PC 1 VS PC 2
g <- ggbiplot(morph_pca, choices = c(1,2), obs.scale = 1, var.scale = 1, 
              groups = morph$SPP, ellipse = T, 
              circle = F,var.axes = T, alpha=0)
{
  g = g +
  #scale_color_manual(name="", values=c("green", "yellow","pink","lightgreen","cyan",
  #                                     "red","darkgreen","magenta","orange","white","purple")) +
  scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
                                       "red","cyan","red","white","white","white")) +
  #scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
  #                                     "cyan","red","cyan","cyan","red","white")) +
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(aes(colour=morph$SPP, shape=morph$SPP))
  g = g + ggtitle("")
  g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'bottom',
                 legend.background = element_rect(fill = 'black', colour = 'white'),
                 legend.text = element_text(colour = 'white'),
                 legend.key = element_rect(fill="black",colour="black"))
  g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
  g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
                panel.grid = element_blank())
  g = g + theme(axis.text.x=element_text(colour="white"),
                axis.text.y=element_text(colour="white"))
  g = g + theme(plot.title = element_text(color="white"),
                axis.title.x = element_text(color="white"),
                axis.title.y = element_text(color="white"))
  print(g)
}

png("~/Documents/Classes/Stephanie/morph_pca_1vs2_spp_arrows_str_26feb2017.png",
    width=12,height=6,pointsize=1,units="cm",
    res=500)
#par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()



## PC 1 VS PC 3
g <- ggbiplot(morph_pca, choices = c(1,3), obs.scale = 1, var.scale = 1, 
              groups = morph$SPP, ellipse = T, 
              circle = F,var.axes = T, alpha=0)
{
  g = g +
  #scale_color_manual(name="", values=c("green", "yellow","pink","lightgreen","cyan",
  #                                     "red","darkgreen","magenta","orange","white","purple")) +
  scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
                                       "red","cyan","red","white","white","white")) +
  #scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
  #                                     "cyan","red","cyan","cyan","red","white")) +
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(aes(colour=morph$SPP, shape=morph$SPP))
  g = g + ggtitle("")
  g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'bottom',
                 legend.background = element_rect(fill = 'black', colour = 'white'),
                 legend.text = element_text(colour = 'white'),
                 legend.key = element_rect(fill="black",colour="black"))
  g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
  g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
                panel.grid = element_blank())
  g = g + theme(axis.text.x=element_text(colour="white"),
                axis.text.y=element_text(colour="white"))
  g = g + theme(plot.title = element_text(color="white"),
                axis.title.x = element_text(color="white"),
                axis.title.y = element_text(color="white"))
  print(g)
}

png("~/Documents/Classes/Stephanie/morph_pca_1vs3_spp_arrows_str_26feb2017.png",
    width=12,height=6,pointsize=1,units="cm",
    res=500)
#par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()


g <- ggbiplot(morph_pca, choices = c(2,3), obs.scale = 1, var.scale = 1, 
              groups = as.character(morph$Structure), ellipse = T, 
              circle = F,var.axes = T, alpha=0)
{
  g = g +
  #scale_color_manual(name="", values=c("green", "yellow","pink","lightgreen","cyan",
  #                                     "red","darkgreen","magenta","orange","white","purple")) +
  #scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
  #                                     "red","cyan","red","white","white","white")) +
  scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
                                       "cyan","red","cyan","cyan","red","white")) +
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(aes(colour=as.character(morph$Structure), shape=as.character(morph$Structure)))
  g = g + ggtitle("")
  g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'bottom',
                 legend.background = element_rect(fill = 'black', colour = 'white'),
                 legend.text = element_text(colour = 'white'),
                 legend.key = element_rect(fill="black",colour="black"))
  g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
  g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
                panel.grid = element_blank())
  g = g + theme(axis.text.x=element_text(colour="white"),
                axis.text.y=element_text(colour="white"))
  g = g + theme(plot.title = element_text(color="white"),
                axis.title.x = element_text(color="white"),
                axis.title.y = element_text(color="white"))
  print(g)
}


g <- ggbiplot(morph_pca, choices = c(2,3), obs.scale = 1, var.scale = 1, 
              groups = morph$SPP, ellipse = T, 
              circle = F,var.axes = T, alpha=0)
{
  g = g +
  #scale_color_manual(name="", values=c("green", "yellow","pink","lightgreen","cyan",
  #                                     "red","darkgreen","magenta","orange","white","purple")) +
  scale_color_manual(name="", values=c("cyan", "white","red","cyan","cyan",
                                       "red","cyan","red","white","white","white")) +
  #scale_color_manual(name="", values=c("cyan", "red","white","cyan","white",
  #                                     "cyan","red","cyan","cyan","red","white")) +
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(aes(colour=morph$SPP, shape=morph$SPP))
  g = g + ggtitle("")
  g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'bottom',
                 legend.background = element_rect(fill = 'black', colour = 'white'),
                 legend.text = element_text(colour = 'white'),
                 legend.key = element_rect(fill="black",colour="black"))
  g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
  g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
                panel.grid = element_blank())
  g = g + theme(axis.text.x=element_text(colour="white"),
                axis.text.y=element_text(colour="white"))
  g = g + theme(plot.title = element_text(color="white"),
                axis.title.x = element_text(color="white"),
                axis.title.y = element_text(color="white"))
  print(g)
}

png("~/Documents/Classes/Stephanie/morph_pca_2vs3_spp_arrows_str_26feb2017.png",
    width=12,height=6,pointsize=1,units="cm",
    res=500)
#par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()


{
g <- ggbiplot(morph_pca, obs.scale = 1, var.scale = 1, 
              groups = morph$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "yellow","pink","lightgreen","cyan",
                                       "red","darkgreen","magenta","orange","white","purple")) +  
  scale_shape_manual(name="", values=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_point(aes(colour=morph$STATE, shape=morph$STATE))
g = g + ggtitle("")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

#[11] "BILL.WIDTH..mm."                 "BILL.HEIGHT..mm."               
#[13] "TARSUS.LENGTH..mm."              "WING.LENGTH.TO.PRIMARIES..mm."  
#[15] "WING.LENGTH.TO.SECONDARIES..mm."

par(mfrow=c(3,2))
boxplot(fusca$BILL.LENGTH..mm.~fusca$WHICH.SIDE.OF.CFB,main="Bill L")

boxplot(fusca$BILL.LENGTH..mm.[fusca$WHICH.SIDE.OF.CFB=="CHIHUAHUAN?"],
        fusca$BILL.LENGTH..mm.[fusca$WHICH.SIDE.OF.CFB=="SONORAN"],
        main="Bill L",
        names=c("C","S"))


boxplot(fusca$BILL.WIDTH..mm.~fusca$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(fusca$BILL.HEIGHT..mm.~fusca$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(fusca$TARSUS.LENGTH..mm.~fusca$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(fusca$WING.LENGTH.TO.PRIMARIES..mm.~fusca$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(fusca$WING.LENGTH.TO.SECONDARIES..mm.~fusca$WHICH.SIDE.OF.CFB,main="Wing S")

## ACTING WEIRD
fusca = fusca[!is.na(fusca$BILL.LENGTH..mm.),]
fusca = fusca[fusca$BILL.LENGTH..mm.!="?",]
fusca = fusca[fusca$WING.LENGTH.TO.SECONDARIES..mm.!="?",]
fusca = fusca[!is.na(fusca$WING.LENGTH.TO.PRIMARIES..mm.),]
fusca_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                          TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.,
                   data=fusca,center=T,scale=T,na.action = "na.omit")
summary(fusca_pca)
fusca_with_pca = merge(fusca,fusca_pca$x,by="row.names")
{
g <- ggbiplot(fusca_pca, obs.scale = 1, var.scale = 1, 
              groups = fusca$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "yellow","magenta","lightgreen","cyan")) +  
  scale_shape_manual(name="", values=c(1,2,3,4,5)) +
  geom_point(aes(colour=fusca$STATE, shape=fusca$STATE))
g = g + ggtitle("Melozone fusca")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/fusca_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
#par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

fusca = fusca[fusca$STATE %in% c("ARIZONA","Arizona","TEXAS","Texas"),]
fusca = fusca[!is.na(fusca$WING.LENGTH.TO.PRIMARIES..mm.),]
fusca_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                     TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.,
                   data=fusca,center=T,scale=T,na.action = "na.omit")
summary(fusca_pca)
fusca_with_pca = merge(fusca,fusca_pca$x,by="row.names")
{
g <- ggbiplot(fusca_pca, obs.scale = 1, var.scale = 1, 
              groups = fusca$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan","magenta","lightgreen","yellow")) +  
  scale_shape_manual(name="", values=c(1,2,3,4,5)) +
  geom_point(aes(colour=fusca$STATE, shape=fusca$STATE))
g = g + ggtitle("Melozone fusca")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/fusca_pca_nomiddle_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
#par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

par(mfrow=c(3,2))
boxplot(cardinalis$BILL.LENGTH..mm.~cardinalis$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(cardinalis$BILL.WIDTH..mm.~cardinalis$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(cardinalis$BILL.HEIGHT..mm.~cardinalis$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(cardinalis$TARSUS.LENGTH..mm.~cardinalis$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(cardinalis$WING.LENGTH.TO.PRIMARIES..mm.~cardinalis$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(cardinalis$WING.LENGTH.TO.SECONDARIES..mm.~cardinalis$WHICH.SIDE.OF.CFB,main="Wing S")

cardinalis_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                          TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                          WING.LENGTH.TO.SECONDARIES..mm.,data=cardinalis,
                        center=T,scale=T)
summary(cardinalis_pca)
cardinalis_with_pca = merge(cardinalis,cardinalis_pca$x,by="row.names")
{
g <- ggbiplot(cardinalis_pca, obs.scale = 1, var.scale = 1, 
              groups = cardinalis$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=cardinalis$STATE, shape=cardinalis$STATE))
g = g + ggtitle("Cardinalis cardinalis")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

boxplot(cardinalis_with_pca$PC1~cardinalis_with_pca$STATE)
cardinalis_with_pca$PC2

png("~/Documents/Classes/Stephanie/cardinalis_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

plot(cardinalis_with_pca$PC1,cardinalis_with_pca$PC2)

par(mfrow=c(3,2))
boxplot(sinuatus$BILL.LENGTH..mm.~sinuatus$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(sinuatus$BILL.WIDTH..mm.~sinuatus$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(sinuatus$BILL.HEIGHT..mm.~sinuatus$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(sinuatus$TARSUS.LENGTH..mm.~sinuatus$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(sinuatus$WING.LENGTH.TO.PRIMARIES..mm.~sinuatus$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(sinuatus$WING.LENGTH.TO.SECONDARIES..mm.~sinuatus$WHICH.SIDE.OF.CFB,main="Wing S")

sinuatus_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                          TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                          WING.LENGTH.TO.SECONDARIES..mm.,data=sinuatus,
                        center=T,scale=T)
summary(sinuatus_pca)
sinuatus_with_pca = merge(sinuatus,sinuatus_pca$x,by="row.names")
{
g <- ggbiplot(sinuatus_pca, obs.scale = 1, var.scale = 1, 
              groups = sinuatus$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=sinuatus$STATE, shape=sinuatus$STATE))
g = g + ggtitle("Cardinalis sinuatus")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}
  
png("~/Documents/Classes/Stephanie/sinuatus_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

par(mfrow=c(3,2))
boxplot(flaviceps$BILL.LENGTH..mm.~flaviceps$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(flaviceps$BILL.WIDTH..mm.~flaviceps$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(flaviceps$BILL.HEIGHT..mm.~flaviceps$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(flaviceps$TARSUS.LENGTH..mm.~flaviceps$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(flaviceps$WING.LENGTH.TO.PRIMARIES..mm.~flaviceps$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(flaviceps$WING.LENGTH.TO.SECONDARIES..mm.~flaviceps$WHICH.SIDE.OF.CFB,main="Wing S")

flaviceps_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                         TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                         WING.LENGTH.TO.SECONDARIES..mm.,data=flaviceps,
                       center=T,scale=T)
summary(flaviceps_pca)
flaviceps_with_pca = merge(flaviceps,flaviceps_pca$x,by="row.names")
{
g <- ggbiplot(flaviceps_pca, obs.scale = 1, var.scale = 1, 
              groups = flaviceps$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=flaviceps$STATE, shape=flaviceps$STATE))
g = g + ggtitle("Auriparus flaviceps")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/flaviceps_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

par(mfrow=c(3,2))
boxplot(brunneicapillus$BILL.LENGTH..mm.~brunneicapillus$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(brunneicapillus$BILL.WIDTH..mm.~brunneicapillus$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(brunneicapillus$BILL.HEIGHT..mm.~brunneicapillus$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(brunneicapillus$TARSUS.LENGTH..mm.~brunneicapillus$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(brunneicapillus$WING.LENGTH.TO.PRIMARIES..mm.~brunneicapillus$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(brunneicapillus$WING.LENGTH.TO.SECONDARIES..mm.~brunneicapillus$WHICH.SIDE.OF.CFB,main="Wing S")

brunneicapillus_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                               TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                               WING.LENGTH.TO.SECONDARIES..mm.,data=brunneicapillus,
                             center=T,scale=T)
summary(brunneicapillus_pca)
brunneicapillus_with_pca = merge(brunneicapillus,brunneicapillus_pca$x,by="row.names")
{
g <- ggbiplot(brunneicapillus_pca, obs.scale = 1, var.scale = 1, 
              groups = brunneicapillus$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=brunneicapillus$STATE, shape=brunneicapillus$STATE))
g = g + ggtitle("Campylorhynchus brunneicapillus")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/brunneicapillus_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()


par(mfrow=c(3,2))
boxplot(curvirostre$BILL.LENGTH..mm.~curvirostre$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(curvirostre$BILL.WIDTH..mm.~curvirostre$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(curvirostre$BILL.HEIGHT..mm.~curvirostre$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(curvirostre$TARSUS.LENGTH..mm.~curvirostre$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(curvirostre$WING.LENGTH.TO.PRIMARIES..mm.~curvirostre$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(curvirostre$WING.LENGTH.TO.SECONDARIES..mm.~curvirostre$WHICH.SIDE.OF.CFB,main="Wing S")

curvirostre_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                           TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                           WING.LENGTH.TO.SECONDARIES..mm.,data=curvirostre,
                         center=T,scale=T)
summary(curvirostre_pca)
curvirostre_with_pca = merge(curvirostre,curvirostre_pca$x,by="row.names")
{
g <- ggbiplot(curvirostre_pca, obs.scale = 1, var.scale = 1, 
              groups = curvirostre$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=curvirostre$STATE, shape=curvirostre$STATE))
g = g + ggtitle("Toxostoma curvirostre")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/curvirostre_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

### melanura 

par(mfrow=c(3,2))
boxplot(melanura$BILL.LENGTH..mm.~melanura$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(melanura$BILL.WIDTH..mm.~melanura$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(melanura$BILL.HEIGHT..mm.~melanura$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(melanura$TARSUS.LENGTH..mm.~melanura$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(melanura$WING.LENGTH.TO.PRIMARIES..mm.~melanura$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(melanura$WING.LENGTH.TO.SECONDARIES..mm.~melanura$WHICH.SIDE.OF.CFB,main="Wing S")

melanura = melanura[!is.na(melanura$BILL.WIDTH..mm.),]
melanura_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                           TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                           WING.LENGTH.TO.SECONDARIES..mm.,data=melanura,
                         center=T,scale=T)
summary(melanura_pca)
melanura_with_pca = merge(melanura,melanura_pca$x,by="row.names")
{
g <- ggbiplot(melanura_pca, obs.scale = 1, var.scale = 1, 
              groups = melanura$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=melanura$STATE, shape=melanura$STATE))
g = g + ggtitle("Polioptila melanura")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/melanura_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()


### bellii 

par(mfrow=c(3,2))
boxplot(bellii$BILL.LENGTH..mm.~bellii$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(bellii$BILL.WIDTH..mm.~bellii$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(bellii$BILL.HEIGHT..mm.~bellii$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(bellii$TARSUS.LENGTH..mm.~bellii$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(bellii$WING.LENGTH.TO.PRIMARIES..mm.~bellii$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(bellii$WING.LENGTH.TO.SECONDARIES..mm.~bellii$WHICH.SIDE.OF.CFB,main="Wing S")

bellii_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                        TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                        WING.LENGTH.TO.SECONDARIES..mm.,data=bellii,
                      center=T,scale=T)
summary(bellii_pca)
bellii_with_pca = merge(bellii,bellii_pca$x,by="row.names")
{
g <- ggbiplot(bellii_pca, obs.scale = 1, var.scale = 1, 
              groups = bellii$STATE, ellipse = T, 
              circle = F,var.axes = F, alpha=0)
g = g +
  scale_color_manual(name="", values=c("green", "cyan")) +  
  scale_shape_manual(name="", values=c(1,2)) +
  geom_point(aes(colour=bellii$STATE, shape=bellii$STATE))
g = g + ggtitle("Vireo bellii")
g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'bottom',
               legend.background = element_rect(fill = 'black', colour = 'white'),
               legend.text = element_text(colour = 'white'),
               legend.key = element_rect(fill="black",colour="black"))
g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
              panel.grid = element_blank())
g = g + theme(axis.text.x=element_text(colour="white"),
              axis.text.y=element_text(colour="white"))
g = g + theme(plot.title = element_text(color="white"),
              axis.title.x = element_text(color="white"),
              axis.title.y = element_text(color="white"))
print(g)
}

png("~/Documents/Classes/Stephanie/bellii_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

## nitens

par(mfrow=c(3,2))
boxplot(nitens$BILL.LENGTH..mm.~nitens$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(nitens$BILL.WIDTH..mm.~nitens$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(nitens$BILL.HEIGHT..mm.~nitens$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(nitens$TARSUS.LENGTH..mm.~nitens$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(nitens$WING.LENGTH.TO.PRIMARIES..mm.~nitens$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(nitens$WING.LENGTH.TO.SECONDARIES..mm.~nitens$WHICH.SIDE.OF.CFB,main="Wing S")

nitens_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                      TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                      WING.LENGTH.TO.SECONDARIES..mm.,data=nitens,
                    center=T,scale=T)
summary(nitens_pca)
nitens_with_pca = merge(nitens,nitens_pca$x,by="row.names")
{
  g <- ggbiplot(nitens_pca, obs.scale = 1, var.scale = 1, 
                groups = nitens$STATE, ellipse = T, 
                circle = F,var.axes = F, alpha=0)
  g = g +
    scale_color_manual(name="", values=c("green", "cyan","white")) +  
    scale_shape_manual(name="", values=c(1,2,3)) +
    geom_point(aes(colour=nitens$STATE, shape=nitens$STATE))
  g = g + ggtitle("Phainopepla nitens")
  g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'bottom',
                 legend.background = element_rect(fill = 'black', colour = 'white'),
                 legend.text = element_text(colour = 'white'),
                 legend.key = element_rect(fill="black",colour="black"))
  g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
  g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
                panel.grid = element_blank())
  g = g + theme(axis.text.x=element_text(colour="white"),
                axis.text.y=element_text(colour="white"))
  g = g + theme(plot.title = element_text(color="white"),
                axis.title.x = element_text(color="white"),
                axis.title.y = element_text(color="white"))
  print(g)
}

png("~/Documents/Classes/Stephanie/nitens_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()

## bilineata

par(mfrow=c(3,2))
boxplot(bilineata$BILL.LENGTH..mm.~bilineata$WHICH.SIDE.OF.CFB,main="Bill L")
boxplot(bilineata$BILL.WIDTH..mm.~bilineata$WHICH.SIDE.OF.CFB,main="Bill W")
boxplot(bilineata$BILL.HEIGHT..mm.~bilineata$WHICH.SIDE.OF.CFB,main="Bill H")
boxplot(bilineata$TARSUS.LENGTH..mm.~bilineata$WHICH.SIDE.OF.CFB,main="Tarsus")
boxplot(bilineata$WING.LENGTH.TO.PRIMARIES..mm.~bilineata$WHICH.SIDE.OF.CFB,main="Wing P")
boxplot(bilineata$WING.LENGTH.TO.SECONDARIES..mm.~bilineata$WHICH.SIDE.OF.CFB,main="Wing S")

bilineata_pca = prcomp(~ BILL.LENGTH..mm.+BILL.WIDTH..mm.+BILL.HEIGHT..mm.+
                      TARSUS.LENGTH..mm.+WING.LENGTH.TO.PRIMARIES..mm.+
                      WING.LENGTH.TO.SECONDARIES..mm.,data=bilineata,
                    center=T,scale=T)
summary(bilineata_pca)
bilineata_with_pca = merge(bilineata,bilineata_pca$x,by="row.names")
{
  g <- ggbiplot(bilineata_pca, obs.scale = 1, var.scale = 1, 
                groups = bilineata$STATE, ellipse = T, 
                circle = F,var.axes = F, alpha=0)
  g = g +
    scale_color_manual(name="", values=c("green", "cyan","white")) +  
    scale_shape_manual(name="", values=c(1,2,3)) +
    geom_point(aes(colour=bilineata$STATE, shape=bilineata$STATE))
  g = g + ggtitle("Amphispiza bilineata")
  g = g + theme(plot.title = element_text(color="white", face="italic",hjust=0))
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'bottom',
                 legend.background = element_rect(fill = 'black', colour = 'white'),
                 legend.text = element_text(colour = 'white'),
                 legend.key = element_rect(fill="black",colour="black"))
  g = g + theme(panel.background = element_rect(fill = 'black', colour = 'grey'))
  g = g + theme(plot.background = element_rect(fill = 'black', colour = 'grey'),
                panel.grid = element_blank())
  g = g + theme(axis.text.x=element_text(colour="white"),
                axis.text.y=element_text(colour="white"))
  g = g + theme(plot.title = element_text(color="white"),
                axis.title.x = element_text(color="white"),
                axis.title.y = element_text(color="white"))
  print(g)
}

png("~/Documents/Classes/Stephanie/bilineata_pca_26feb2017.png",
    width=12,height=8,pointsize=1,units="cm",
    res=500)
par(mfrow=c(1,1),mar=c(0,0,0,0))
print(g)
dev.off()


###

# mean and variation of variables by species 

species = as.character(unique(morph_with_pca$SPP))
variables = colnames(morph_with_pca)[c(11:16,28:34)]
varnum = c(11:16,28:34)

for(s in species){
  for(v in variables){
    d = morph_with_pca[morph_with_pca$SPP==s,]
    summ = (summary(d[,v]))
    print(c(s,v,summ))
  }
}

##
mean(morph_with_pca$PC1[morph_with_pca$SPP=="flaviceps"])
mean(morph_with_pca$PC1[morph_with_pca$SPP=="bellii"])
mean(morph_with_pca$PC1[morph_with_pca$SPP=="curvirostre"])
mean(morph_with_pca$PC1[morph_with_pca$SPP=="cardinalis"])
mean(morph_with_pca$PC1[morph_with_pca$SPP=="melanura"])
mean(morph_with_pca$PC1[morph_with_pca$SPP=="brunneicapillus"])
mean(morph_with_pca$PC1[morph_with_pca$SPP=="fusca"])

spp = c("flaviceps","bellii","curvirostre","cardinalis","melanura","brunneicapillus","fusca","bilineata","sinuatus","nitens")
means = c(mean(morph_with_pca$PC1[morph_with_pca$SPP=="flaviceps"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="bellii"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="curvirostre"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="cardinalis"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="melanura"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="brunneicapillus"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="fusca"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="bilineata"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="sinuatus"]),
          mean(morph_with_pca$PC1[morph_with_pca$SPP=="nitens"])
)

spp = c("flaviceps","bellii","curvirostre","cardinalis","melanura","brunneicapillus","fusca","bilineata","sinuatus","nitens")


means2 = c(mean(morph_with_pca$PC2[morph_with_pca$SPP=="flaviceps"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="bellii"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="curvirostre"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="cardinalis"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="melanura"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="brunneicapillus"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="fusca"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="bilineata"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="sinuatus"]),
           mean(morph_with_pca$PC2[morph_with_pca$SPP=="nitens"])
)


means3 = c(mean(morph_with_pca$PC3[morph_with_pca$SPP=="flaviceps"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="bellii"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="curvirostre"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="cardinalis"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="melanura"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="brunneicapillus"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="fusca"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="bilineata"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="sinuatus"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="nitens"])
)

means3 = c(mean(morph_with_pca$PC3[morph_with_pca$SPP=="flaviceps"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="bellii"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="curvirostre"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="cardinalis"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="melanura"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="brunneicapillus"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="fusca"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="bilineata"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="sinuatus"]),
           mean(morph_with_pca$PC3[morph_with_pca$SPP=="nitens"])
)




"              Pi_Mean       Pi_Var
verdin_pi 0.002412744 3.680729e-06
wren_pi   0.002423780 4.219876e-06
gnat_pi   0.002732944 2.675438e-06
fusca_pi  0.005099060 1.171104e-05
vireo_pi  0.013160039 4.417989e-05
thrash_pi 0.014618547 6.254936e-05
card_pi   0.016618799 6.955303e-05"

pi = c(0.002412744,0.013160039,0.014618547,0.016618799,0.002732944,0.002423780,0.005099060,NA,NA,NA)

"bil_prop20 #0.03841932
bel_prop20 #0.2755214
mel_prop20 #0.2184413
sin_prop20 #0.07025247
fus_prop20 #0.004390779
fla_prop20 #0.06586169
nit_prop20 #0.07903403
bru_prop20 #0.04390779
cur_prop20 #0.03512623"

#suit20 = c("0.06586169","0.2755214","0.03512623","cardinalis","0.2184413","0.04390779","0.004390779","0.03841932","0.07025247","0.07903403")
suit20 = c("0.06586169","0.2755214","0.03512623",NA,"0.2184413","0.04390779","0.004390779","0.03841932","0.07025247","0.07903403")

# amphispiza	0.19570312	
# bellii	0.14259912	
# brunneicapillus   0.12123203	
# curvirostre 0.09289573	
# flaviceps	0.01578953	
# fusca	0.09119378	
# phainopepla	0.18654452	
# polioptila	0.02024653	
# sinuatus	0.05652053
suit = c()
thresh = c("0.01578953","0.14259912","0.09289573",NA,"0.02024653",
        "0.12123203","0.09119378","0.19570312","0.05652053","0.18654452")


"bil_propThresh #0.03841932
bel_propThresh #0.219539
mel_propThresh #0.02744237
sin_propThresh #0.04500549
fus_propThresh #0
fla_propThresh #0.008781559
nit_propThresh #0.07903403
bru_propThresh #0.03732162
cur_propThresh #0.03293085"
suitThresh = c("0.008781559","0.219539","0.03293085",NA,"0.02744237","0.03732162","0","0.03841932","0.04500549","0.07903403")

png("~/Documents/Classes/Stephanie/pc1_vs_pi_26feb2017.png",
    width=12,height=8,pointsize=12,units="cm",
    res=500)
par(mfrow=c(1,1),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
plot(pi,means,col=c("white"),
     pch = c(1),xlim=c(0,0.017),ylim=c(-3,3),
     ylab="Mean PC1",xlab="Nucleotide Diversity")
axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
abline(lm(means~pi),col="white")
summary(lm(means~pi))
dev.off()

plot(pi,means2,col=c("white"),
     pch = c(1),xlim=c(0,0.017),ylim=c(-2,3),
     ylab="Mean PC2",xlab="Nucleotide Diversity")
axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
abline(lm(means2~pi),col="white")
summary(lm(means2~pi))

plot(pi,means3,col=c("white"),
     pch = c(1),xlim=c(0,0.017),ylim=c(-0.5,0.5),
     ylab="Mean PC3",xlab="Nucleotide Diversity")
axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
abline(lm(means3~pi),col="white")
summary(lm(means3~pi))

png("~/Documents/Classes/Stephanie/suitability20_vs_pi_26feb2017.png",
    width=12,height=8,pointsize=12,units="cm",
    res=500)
par(mfrow=c(1,1),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
plot(pi,suit20,col=c("white"),
     pch = c(1),xlim=c(0,0.017),ylim=c(0,0.3),
     ylab="Proportion Unsuitable",xlab="Nucleotide Diversity")
axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
#abline(lm(suit20~pi),col="white")
summary(lm(suit20~pi))
dev.off()

png("~/Documents/Classes/Stephanie/suitThresh_vs_pi_26feb2017.png",
    width=12,height=8,pointsize=12,units="cm",
    res=500)
par(mfrow=c(1,1),bg="black",col.axis="white",col.lab="white",
    col.main="white",col.sub="white",
    lwd=1,font.main=4,mar=c(2,2,2,1))
plot(pi,suitThresh,col=c("white"),
     pch = c(1),xlim=c(0,0.017),
     ylab="Proportion Unsuitable",xlab="Nucleotide Diversity")
axis(1,col="white",col.ticks="white",lwd=2); axis(2,col="white",col.ticks="white",lwd=2)
#abline(lm(suitThresh~pi),col="white")
summary(lm(suitThresh~pi))
dev.off()
