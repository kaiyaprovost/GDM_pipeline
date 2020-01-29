## TO DO: go through once with tutorial data

## dismo setup for mac
Sys.setenv(NOAWT=TRUE) #only for certain macs before loading rJava
library(rJava) # required to run Maxent within dismo
options(java.parameters = "-Xmx1g" ) #optional, to set the memory allocated to java, hence maxent (has to be done before loading dismo)
library(dismo) # should also load automatically the required packages sp and raster
library(rgdal)
library(raster)

## crop worldclim?
#cropfiles = list.files("/Users/kprovost/Documents/Classes/GIS/Worldclim/rawWorldclimTif",pattern=".tif$",
#                       full.names=FALSE)
#e<-extent(c(-167, -55, 6, 73))
# crop bioclim data to this extent
#for(file in cropfiles){
#  test = raster(paste("/Users/kprovost/Documents/Classes/GIS/Worldclim/rawWorldclimTif/",file,sep=""))
#  croptest <-crop(test, e)
#  writeRaster(croptest,filename=paste("/Users/kprovost/Documents/Classes/GIS/Worldclim/crop_",file,".asc",sep=""),
#              overwrite=FALSE,format="ascii")
#}

## get and stack rasters
filenames = list.files("/Users/kprovost/Documents/Classes/GIS/Worldclim",pattern="crop*",
                       full.names=TRUE)
## TIF would be better!
#filenames
#worldclim = stack(filenames)
#worldclim
#writeRaster(worldclim, filename="/Users/kprovost/Documents/Classes/GIS/Worldclim/crop_bioclim_all_stack.tif", options="INTERLEAVE=BAND", 
#            overwrite=TRUE)

#pred_nf = worldclim
pred_nf = stack("/Users/kprovost/Documents/Classes/GIS/Worldclim/crop_bioclim_all_stack.tif")

runNicheModels = function(path,csvname,pred_nf,filenames,outfilename){
  ## import thinned datasets 
  print("import")
  vireo = read.csv(paste(path,csvname,sep=""))
  vireo = vireo[,c(2,3)]
  names(vireo) = c("longitude","latitude")
  
  ## visualize setophaga
  library(maptools)
  data(wrld_simpl)
  #plot(wrld_simpl,xlim=c(-113,-96),ylim=c(20,34),
  #     axes=TRUE,col="light yellow")
  #points(vireo$longitude,vireo$latitude,
  #       col="orange",pch=20,cex=0.75,
  #       xlim=c(-113,-96),ylim=c(20,34))
  
  ## create random bg points
  print("bg")
  mask = raster(filenames[1])
  set.seed(120412)
  background = randomPoints(mask,500)
  colnames(background)= c("latitude","longitude")
 
  ## extract seto values and bg values
  print("presvals")
  presvals = extract(pred_nf, vireo)
  head(presvals)
  bgvals = extract(pred_nf,background)
  pb = c(rep(1,nrow(presvals)),rep(0,nrow(bgvals)))
  sdmdata = data.frame(cbind(pb,rbind(presvals,bgvals)))
  
  ## do a kfold
  print("kfold")
  k = 5
  group = kfold(vireo,5)
  group[1:10]
  unique(group)
  e = list()
  for(i in 1:k){
    train = vireo[group != i,]
    print("1")
    test = vireo[group == i,]
    print("2")
    bc = bioclim(train)
    print("3")
    e[[i]] = evaluate(p=test,a=background,bc) 
  }
  e
  auc = sapply(e,function(x){slot(x,"auc")})
  tpr = sapply(e,function(x){x@t[which.max(x@TPR+x@TNR)]})
  print(auc)
  print(tpr)
  
  ## do the kfold, but only once
  print("kfold 1")
  group = kfold(vireo,5)
  pres_train = vireo[group != 1,]
  pres_test = vireo[group == 1,]
  group = kfold(background,5)
  back_train = background[group!=1,]
  back_test = background[group==1,]
  
  r = raster(pred_nf,1)
  png(filename=paste(path,outfilename,"_traintest.png",sep="")) 
  {
  plot(!is.na(r),col=c("white","light grey"),legend=FALSE)
  points(back_train, pch='-', cex=0.5, col='yellow')
  points(back_test, pch='-', cex=0.5, col='black')
  points(pres_train, pch= '+', col='green')
  points(pres_test, pch='+', col='blue')
  } 
  dev.off()
  
  ## BIOCLIM MODEL
  print("bioclim")
  bc = bioclim(pred_nf,pres_train)
  #plot(bc,a=1,b=2,p=0.85)
  e = evaluate(pres_test,back_test,bc,pred_nf)
  e
  tr = threshold(e,"spec_sens")
  tr
  pb = predict(pred_nf,bc,progress="")
  pb
  
  writeRaster(pb,filename=paste(path,outfilename,"_threshval_",as.character(tr),"_bioclim.asc",sep=""),overwrite=T)
  
  
  png(filename=paste(path,outfilename,"_bioclim.png",sep="")) 
  {
  par(mfrow=c(1,2))
  plot(pb, main='Bioclim, raw values')
  plot(wrld_simpl, add=TRUE, border='dark grey')
  plot(pb > tr, main='presence/absence')
  plot(wrld_simpl, add=TRUE, border='dark grey')
  points(pres_train, pch='+')
  } 
  dev.off()
  
  ## MAXENT MODEL
  print("maxent")
  system.file("java", package="dismo")
  
  ## NEED TO INSTALL MAXE
  
  jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
  if (file.exists(jar)) {
    xm <- maxent(pred_nf, pres_train)
    print("done xm")
    #plot(xm)
    response(xm)
    e = evaluate(pres_test,back_test,xm,pred_nf)
    e
    px = predict(pred_nf,xm,progress="")
    png(filename=paste(path,outfilename,"_maxent.png",sep=""))
    {
    par(mfrow = c(1,2))
    plot(px,main="Maxent, raw values")
    plot(wrld_simpl,add=TRUE,border="dark grey")
    tr = threshold(e,"spec_sens")
    plot(px > tr, main = "pres/absence")
    plot(wrld_simpl,add=TRUE,border="dark grey")
    points(pres_train,pch="+")
    }
    dev.off()
    
    writeRaster(px,filename=paste(path,outfilename,"_maxent_threshval_",tr,"_.asc",sep=""),overwrite=T)
    
    
  } else {
    print(cat('cannot run this example because maxent is not available'))
    #plot(1)
  }
   
}

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Polioptila melanuraAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_18_09/",
#                          csvname="Polioptila melanuraAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_18_09_thin1.csv",
#                          pred_nf=pred_nf,
#                          filenames=filenames,
#                          outfilename="polioptila_tr")

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Vireo belliiAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26/",
#                          csvname="Vireo belliiAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26_thin1.csv",
#                          pred_nf=pred_nf,
#                          filenames=filenames,
#                          outfilename="bellii_tr")

##

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Phainopepla nitensAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_43_40/",
#               csvname="Phainopepla nitensAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_43_40_thin1.csv",
#               pred_nf=pred_nf,
#               filenames=filenames,
#               outfilename="phainopepla_tr")

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Melozone fuscaAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_04_38_29/",
#               csvname="Melozone fuscaAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_04_38_29_thin1.csv",
#               pred_nf=pred_nf,
#               filenames=filenames,
#               outfilename="fusca_tr")

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Cardinalis sinuatusAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_46_31/",
#               csvname="Cardinalis sinuatusAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_27_20_46_31_thin1.csv",
#               pred_nf=pred_nf,
#               filenames=filenames,
#               outfilename="sinuatus_tr")

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Auriparus flavicepsAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26/",
#               csvname="Auriparus flavicepsAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26_thin1.csv",
#               pred_nf=pred_nf,
#               filenames=filenames,
#               outfilename="flaviceps_tr")

##

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Amphispiza bilineataAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_00_16_23/",
#               csvname="Amphispiza bilineataAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_00_16_23_thin1.csv",
#               pred_nf=pred_nf,
#               filenames=filenames,
#               outfilename="amphispiza_tr")

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Campylorhynchus brunneicapillusAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26/",
#               csvname="Campylorhynchus brunneicapillusAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_06_42_26_thin1.csv",
#               pred_nf=pred_nf,filenames=filenames,
#               outfilename="brunneicapillus")

#runNicheModels(path="/Users/kprovost/Documents/Dissertation/GBIF/Toxostoma curvirostreAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_03_11_12/",
#               csvname="Toxostoma curvirostreAllSpecies_NotThinned_NoUncertains_WITHINPOLYGON_out_2017_10_28_03_11_12_thin1.csv",
#               pred_nf=pred_nf,filenames=filenames,
#               outfilename="curvirostre")


