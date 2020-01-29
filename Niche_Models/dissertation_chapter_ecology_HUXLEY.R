## CHANGE THIS TO DETERMINE WHERE TO START
FIRSTSTEP = 7

#detach("package:subsppLabelR", unload = TRUE)
library(devtools)
devtools::install_github('kaiyaprovost/subsppLabelR', force = F)
devtools::install_github("cran/ecospat", force = F)
devtools::install_github("danlwarren/ENMTools", force = F)
devtools::install_github("rsh249/rasterExtras", force = F)
devtools::install_github("petrelharp/landsim", force = F)
library(subsppLabelR)
library(ecospat)
library(ENMTools)
library(rasterExtras)
library(landsim)
library(sys)

#sudo R CMD javareconf JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/A/Headers/

dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package, ")"
  ))))
}

packages = c(
  "dismo",
  "GISTools",
  "raster",
  "rgdal",
  "ENMeval",
  "phyloclim",
  "data.table",
  "dplyr",
  "EMCluster",
  "knor",
  "maps",
  "MASS",
  "parallel",
  "plotly",
  "rgeos",
  "roxygen2",
  "rworldmap",
  "sp",
  "spThin",
  "spocc",
  "spThin",
  "viridis",
  "auk",
  "rebird"
)
## rJava is causing issues

for (p in packages) {
  dynamic_require(p)
}

if (FIRSTSTEP <= 7) {
  
  plotStuff=T
  
  Env = stack("/home/kprovost/nas1/ENMs/ENMS_multilayer4.tif")
  
  thinDirectory = "/home/kprovost/nas1/ENMs/THIN/"
  setwd(thinDirectory)
  
  thinned_csvs = list.files(
    path = thinDirectory,
    pattern = "thin1_USA.csv$",
    full.names = T
  )

  outdirectory="/home/kprovost/nas1/ENMs/USA_ONLY/"
  setwd(outdirectory)
  
  for (thinfile in thinned_csvs) {
    
    print(thinfile)
    thinned = read.csv(thinfile)
    sppname = thinned$name[1]
    loc = cbind(thinned$longitude,thinned$latitude)
    extr = extract(Env, loc) ## get env data at each location
    loc = loc[!is.na(extr[,1]),] ## remove any locations that have no ENV data
    
    ## add in more of this:
    ## https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html
    
    ## split data into four equal blocks and different feature classes 
    ## throws errors with hinge 
    res = ENMevaluate(occ=loc, env = Env, method='block', 
                      parallel=T, numCores=8, fc=c("L", "LQ"), 
                      RMvalues=seq(0.5,4,0.5), rasterPreds=T,
                      updateProgress = T)
    ## this uses your occurences and trimmed env data. it will
    ## separate the data into 4 blocks, split sample into 4 equal quads
    ## running in parallel, with 4 cores. 
    ## the feature classes used are Linear, Linear+Quadratic, and Hunge.
    ## the regularization values will range from 1/2 to 4 by 1/2s.
    ## and rasterPreds=T means not running the predict function we ran before
    ## and thus can't have any AICc
    
    ## now find the best models 
    ## minimize omission rate and optimize AUC values
    ## remeber to change ORmin to OR10 depending on sample size
    ## for omission rate, if small sample size and high confidence
    ## do mean.ORmin. if don't know, do mean.OR10
    setsort = res@results[order(res@results[,'avg.test.or10pct']),] ## previously Mean.ORmin
    setsort2 = setsort[order(setsort[,'avg.test.AUC'], decreasing=TRUE),] ## previously Mean.AUC
    top = setsort2[1,]
    print(top)
    write.csv(setsort2,file=paste("AddedPredThin_ResultsTable_",sppname,"_addedLayers_zoomedin.csv",sep=""))
    
    best = which(as.character(res@results[,1]) == as.character(setsort2[1,1]))
    
    if(plotStuff==T) {
      pred.raw = predict(Env, res@models[[best]])
      writeRaster(pred.raw,filename=paste("AddedPredThin_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      png(paste("PredThin_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(pred.raw, col=viridis::viridis(99),main=sppname)
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      
      png(paste("PredThin_avg.test.AUC_",sppname,"_addedLayers_zoomedin.png",sep=""))
      eval.plot(res@results, "avg.test.AUC")
      dev.off()
    }
    
    ## with thresh model
    ev.set <- evaluate(thinned[,2:3], res@bg.pts, res@models[[best]], Env)
    th1 = threshold(ev.set) ## omission options
    
    write.csv(th1,file=paste("AddedPredThin_ThreshTable_",sppname,"_addedLayers_zoomedin.csv",sep=""))
    
    if(plotStuff==T) {
      
      p1.nomit = pred.raw>= th1$no_omission ## the highest thresh where no points are omittec
      p1.equal = pred.raw>= th1$equal_sens_spec ## equal sensitivity and specificity according to ROC curve
      p1.spsen = pred.raw>= th1$spec_sens ## equal sensitivity and specificity according to ROC curve
      
      writeRaster(p1.nomit,filename=paste("Thresh_NoOmission_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      writeRaster(p1.equal,filename=paste("Thresh_EqualSensSpec_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      writeRaster(p1.spsen,filename=paste("Thresh_MaxSpecSens_BestModel_",sppname,"_addedLayers_zoomedin.asc",sep=""),
                  format="ascii",overwrite=T)
      png(paste("Thresh_NoOmission_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p1.nomit, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      png(paste("Thresh_EqualSensSpec_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p1.equal, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
      png(paste("Thresh_MaxSpecSens_BestModel_",sppname,"_addedLayers_zoomedin.png",sep=""))
      plot(p1.spsen, col=viridis::viridis(99))
      points(loc,pch=3,col=rgb(1,0,0,0.1))
      dev.off()
    }
    
    ## this is where it was getting stuck 
    ## remove stuff from memory now that done 
    
    rm(p1.equal)
    rm(p1.nomit)
    rm(p1.spsen)
    rm(pred.raw)
    rm(res)
    rm(setsort)
    rm(setsort2)
    rm(th1)
    rm(thinned)
    rm(top)
    rm(loc)
    rm(extr)
    rm(best)
    rm(ev.set)
    
  }
  
}

