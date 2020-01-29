library(scales)

for (p in packages) {
  dynamic_require(p)
}

#path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/USA_ONLY/"
#path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WITH_MEXICO/"
path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/"

setwd(path)

#files = list.files(pattern="AddedPredThin.+USA.asc$",recursive = T)
#files = list.files(pattern="AddedPredThin.+zoomedin.asc$",recursive = T)
#files = list.files(pattern="AddedPredThin.+addedLAyers.asc$",recursive = T)
#files = list.files(pattern="AddedPredThin.+LGM.asc$",recursive = T)
#files = list.files(pattern="AddedPredThin.+MID.asc$",recursive = T)
files = list.files(pattern="AddedPredThin.+worldclim.asc$",recursive = T)

doThresh=T

for (asc in files) {
  print(asc)
  print("--reading")
  ras = raster(asc)
  print("--calculating")
  maximum = max(values(ras),na.rm=T)
  minimum = min(values(ras),na.rm=T)
  print("--scaling raster")
  values(ras) = scales::rescale(values(ras),to=c(0,1))
  newasc = paste(asc,"_",maximum,"_",minimum,".rescaled",sep="")
  print("--writing raster")
  writeRaster(ras,newasc,format="ascii",overwrite=T)
  
  if(doThresh==T) {
  
  print("--reading thresh") 
  tab = sub(pattern="BestModel",replacement="ThreshTable",x=asc)
  tab = sub(pattern="asc",replacement="csv",x=tab)
  newtab = paste(tab,".rescaled.csv",sep="")
  read = read.csv(tab)
  print("--scaling thresh")
  newread = c("threshold",as.vector(sapply(read[,2:7],FUN=function(x){(x-minimum)/(maximum-minimum)})))
  names(newread) = names(read)
  names(newread)[1] = ""
  newread = rbind(newread)
  print("--writing thresh")
  write.csv(newread,newtab,row.names=F)
  }
}
