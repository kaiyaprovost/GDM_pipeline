#devtools::install_github("wpeterman/ResistanceGA", build_vignettes = F,force=F) # Download package
library(ResistanceGA)
library(GA)
library(gdistance)
library(scales)
rm(list = ls())

# Rscript "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/resistanceGA_stuff.R"

aggregationFactor = 1
dobigraster = F
dosmallraster = T
doabund = T
dosmallnospp = T

outputdirectory = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/"
#outputdirectory = "/home/kprovost/nas1/ENMs/DISTANCES/"
setwd(outputdirectory)

#samples = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_BY_SPECIES.txt",sep="\t")
#samples = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_ONLY.txt",sep="\t")

#samples = read.csv("/home/kprovost/nas1/ENMs/LATLONGS_BY_SPECIES.txt",sep="\t")
#samples = read.csv("/home/kprovost/nas1/ENMs/LATLONGS_ONLY.txt",sep="\t")

#allspp = as.character(sort(unique(samples$spp)))
#print(allspp)

## now do the same thing but for the entire dataset at once 
## with full lgm etc stuff 

if (dobigraster == T) {
  
  samples = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_BY_SPECIES.txt",sep="\t")
  allspp = as.character(sort(unique(samples$spp)))
  allspp = allspp[allspp!="CARDINALIS"]
  
  whereraster="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/"
  #whereraster="/home/kprovost/nas1/ENMs/ORIGINALS/"
  
  pathlist = list.files(path=whereraster,
                        pattern=glob2rx(paste("ENMS","*.tif",sep="")),
                        #pattern=glob2rx(paste("*.asc",sep="")),
                        recursive = T)
  #pathlist = pathlist[lapply(pathlist,function(x) length(grep("USA",x,value=FALSE))) == 0]
  #pathlist = pathlist[lapply(pathlist,function(x) length(grep("zoomed",x,value=FALSE))) == 0]
  
  sample.locales = (samples[,c("long","lat")])
  #sample.locales = (samples[,c(3,2)])
  #sample.locales[,1] = as.numeric(levels(sample.locales[,1]))[sample.locales[,1]]
  #sample.locales[,2] = as.numeric(levels(sample.locales[,2]))[sample.locales[,2]]
  sample.locales = unique(sample.locales)
  sample.locales = sample.locales[complete.cases(sample.locales),]
  
  maxlat = max(sample.locales$lat,na.rm = T)+1
  maxlong = max(sample.locales$long,na.rm = T)+1
  minlat = min(sample.locales$lat,na.rm = T)-1
  minlong = min(sample.locales$long,na.rm = T)-1
  
  sample.locales = unique(sample.locales)
  sample.locales = SpatialPoints(sample.locales)
  
  
  print("setting extent")
  ext = extent(minlong,maxlong,minlat,maxlat)
  print(ext)
  
  for (i in 1:length(pathlist)) {
    
    print(paste(i,"of",length(pathlist)))
    old <- Sys.time()
    print(old)
    
    rasterpath = paste(whereraster,pathlist[i],sep="")
    #rasterpath = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc_-6.46640157699585_-69.5545196533203.rescaled.asc"
    
    print(rasterpath)
    
    splitpath = strsplit(rasterpath,"/")[[1]]
    suffix = splitpath[length(splitpath)]
    
    #print("load raster")
    continuous = stack(rasterpath)
    numlayers = length(continuous@layers)
    
    
    ## crop raster
    print("cropping")
    continuous = crop(continuous,ext)
    png(paste(outputdirectory,suffix,".AGGREGATE-",aggregationFactor,".png",sep=""))
    par(mfrow=c(1,2))
    plot(continuous[[1]])
    
    ## need to rescale the individual layers 
    print("rescaling")
    
    for (k in 1:numlayers) {
      #print(k)
      todo = continuous[[k]]
      vals = values(todo)
      vals2 = rescale(vals)
      values(todo) = vals2
      continuous[[k]] = todo
      
    }
    continuous = stack(continuous)
    
    if (aggregationFactor > 1) {
      
      print(paste("aggregating by factor of",aggregationFactor))
      continuous <- aggregate(continuous, fact=aggregationFactor, fun=mean)
    } else {
      print("skipping aggregation")
      
    }
    
    plot(continuous[[1]])
    points(sample.locales)
    dev.off()
    
    
    #print("making temps")
    temp = sample.locales
    temp = SpatialPoints(temp)
    #temp = temp[!is.na(extract(continuous,temp)),]
    temp = temp[complete.cases(extract(continuous,temp)),]
    
    individuals = as.numeric(rownames(temp@coords))
    included = samples[individuals,]
    #catalog = included$CATALOG.NUMBER
    latlong = paste(included$lat,included$long)
    
    ## get euclidian distances
    
    rawvalues = extract(continuous,temp)
    rownames(rawvalues) = latlong
    rawdist = as.data.frame(as.matrix(dist(rawvalues)))
    rawfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".RAW-EUCLIDIAN.txt",sep="")
    write.csv(rawdist,rawfile)
    
    
    
    ## may need to do GA prep or gdist prep on the stack? 
    
    continuous = stack(continuous)
    
    GA.inputs <- GA.prep(ASCII.dir = continuous,
                         Results.dir = outputdirectory)
    
    outfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".MORPH-AND-GENE-DISTANCES.txt",sep="")
    #print(paste("outfile:",outfile))
    #print(outfile)
    
    #print("gdist prep")
    gdist.inputs <- gdist.prep(length(temp),
                               samples = temp,
                               ##longlat = F, ## may need to comment
                               method = 'costDistance') ## works with costdist?
    
    print("gdist response")
    
    PARM=rep(c(1,3.5,150),times=numlayers)
    
    Resist <- Combine_Surfaces(PARM = PARM,
                               gdist.inputs = gdist.inputs,
                               GA.inputs = GA.inputs,
                               out = NULL,
                               rescale = TRUE)
    values(Resist) = rescale(values(Resist))
    
    png(paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".FULL_RESISTANCE_SURFACE.png",sep=""))
    plot(Resist)
    dev.off()
    
    gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                    #r = continuous,
                                    r=Resist,
                                    scl=F) # works if latlongF and methodCost
    
    #print("matrix")
    distdf = as.matrix(gdist.response)
    #print("rename")
    rownames(distdf) = latlong
    colnames(distdf) = latlong
    
    #print("output")
    write.csv(distdf,outfile)
    new <- Sys.time() - old # calculate difference
    print(paste("Time elapsed:",new))
    
  }
}

if (dosmallraster == T) {
  
  samples = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_BY_SPECIES.txt",sep="\t")
  allspp = as.character(sort(unique(samples$spp)))
  allspp = allspp[allspp!="CARDINALIS"]
  
  whereraster="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/"
  #whereraster="/home/kprovost/nas1/ENMs/RESCALED/"
  
  #for ( j in 6:8) {
  for (j in 1:length(allspp)) {
    print(j)
    
    spp = allspp[j]
    print(spp)
    
    pathlist = list.files(path=whereraster,
                          pattern=glob2rx(paste("*",tolower(spp),"*rescaled*.asc",sep="")),
                          recursive = T)
    pathlist = pathlist[lapply(pathlist,function(x) length(grep("USA",x,value=FALSE))) == 0]
    pathlist = pathlist[lapply(pathlist,function(x) length(grep("zoomed",x,value=FALSE))) == 0]
    
    
    sample.locales = (samples[samples$spp==spp,c(3,2)])
    #sample.locales[,1] = as.numeric(levels(sample.locales[,1]))[sample.locales[,1]]
    #sample.locales[,2] = as.numeric(levels(sample.locales[,2]))[sample.locales[,2]]
    sample.locales = unique(sample.locales)
    sample.locales = sample.locales[complete.cases(sample.locales),]
    
    maxlat = max(sample.locales$lat,na.rm = T)+1
    maxlong = max(sample.locales$long,na.rm = T)+1
    minlat = min(sample.locales$lat,na.rm = T)-1
    minlong = min(sample.locales$long,na.rm = T)-1
    
    print(paste(minlong,maxlong,minlat,maxlat))
    
    print("setting extent")
    ext = extent(minlong,maxlong,minlat,maxlat)
    print(ext)
    
    sample.locales = SpatialPoints(sample.locales)
    
    for (i in 1:length(pathlist)) {
      
      print(paste(i,"of",length(pathlist)))
      old <- Sys.time()
      print(old)
      
      rasterpath = paste(whereraster,pathlist[i],sep="")
      #rasterpath = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc_-6.46640157699585_-69.5545196533203.rescaled.asc"
      
      print(rasterpath)
      
      splitpath = strsplit(rasterpath,"/")[[1]]
      suffix = splitpath[length(splitpath)]
      
      #print("load raster")
      continuous = raster(rasterpath)
      
      ## crop raster
      print("cropping")
      continuous = crop(continuous,ext)
      png(paste(outputdirectory,suffix,".AGGREGATE-",aggregationFactor,".png",sep=""))
      par(mfrow=c(1,2))
      plot(continuous)
      
      
      if (aggregationFactor > 1) {
        
        print(paste("aggregating by factor of",aggregationFactor))
        continuous <- aggregate(continuous, fact=aggregationFactor, fun=mean)
      } else {
        print("skipping aggregation")
        
      }
      
      plot(continuous)
      points(sample.locales)
      dev.off()
      
      temp = sample.locales
      temp = SpatialPoints(temp)
      temp = temp[complete.cases(extract(continuous,temp)),]
      
      individuals = as.numeric(rownames(temp@coords))
      included = samples[individuals,]
      #catalog = included$CATALOG.NUMBER
      latlong = paste(included$lat,included$long)
      
      rawvalues = as.data.frame(extract(continuous,temp))
      rownames(rawvalues) = latlong
      rawdist = as.data.frame(as.matrix(dist(rawvalues)))
      rawfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".RAW-EUCLIDIAN.txt",sep="")
      write.csv(rawdist,rawfile)
      
      
      outfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".MORPH-AND-GENE-DISTANCES.txt",sep="")
      #print(paste("outfile:",outfile))
      #print(outfile)
      
      #print("gdist prep")
      gdist.inputs <- gdist.prep(length(temp),
                                 samples = temp,
                                 ##longlat = F, ## may need to comment
                                 method = 'costDistance') ## works with costdist?
      
      print("gdist response")
      gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                      r = continuous,
                                      scl=F) # works if latlongF and methodCost
      
      #print("matrix")
      distdf = as.matrix(gdist.response)
      #print("rename")
      rownames(distdf) = latlong
      colnames(distdf) = latlong
      
      #print("output")
      write.csv(distdf,outfile)
      new <- Sys.time() - old # calculate difference
      print(paste("Time elapsed:",new))
      
    }
  }
}

if (doabund == T) {
  
  samples = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_ONLY.txt",sep="\t")
  
  #whereraster="/home/kprovost/nas1/ENMs/ABUNDANCE/"
  whereraster="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/abundance_asciis/"
  
  pathlist = list.files(path=whereraster,
                        pattern=glob2rx(paste("*.asc",sep="")),
                        recursive = T)
  
  sample.locales = (samples[,c("long","lat")])
  #sample.locales = (samples[,c(3,2)])
  #sample.locales[,1] = as.numeric(levels(sample.locales[,1]))[sample.locales[,1]]
  #sample.locales[,2] = as.numeric(levels(sample.locales[,2]))[sample.locales[,2]]
  sample.locales = unique(sample.locales)
  sample.locales = sample.locales[complete.cases(sample.locales),]
  
  maxlat = max(sample.locales$lat,na.rm = T)+1
  maxlong = max(sample.locales$long,na.rm = T)+1
  minlat = min(sample.locales$lat,na.rm = T)-1
  minlong = min(sample.locales$long,na.rm = T)-1
  
  sample.locales = unique(sample.locales)
  sample.locales = SpatialPoints(sample.locales)
  
  
  print("setting extent")
  ext = extent(minlong,maxlong,minlat,maxlat)
  print(ext)
  
  for (i in 1:length(pathlist)) {
    
    print(paste(i,"of",length(pathlist)))
    old <- Sys.time()
    print(old)
    
    rasterpath = paste(whereraster,pathlist[i],sep="")
    #rasterpath = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc_-6.46640157699585_-69.5545196533203.rescaled.asc"
    
    print(rasterpath)
    
    splitpath = strsplit(rasterpath,"/")[[1]]
    suffix = splitpath[length(splitpath)]
    
    #print("load raster")
    continuous = raster(rasterpath)
    
    ## crop raster
    print("cropping")
    continuous = crop(continuous,ext)
    png(paste(outputdirectory,suffix,".AGGREGATE-",aggregationFactor,".png",sep=""))
    par(mfrow=c(1,2))
    plot(continuous[[1]])
    
    ## need to rescale the individual layers 
    print("rescaling")
    
    todo = continuous
    vals = values(todo)
    vals2 = rescale(vals)
    values(todo) = vals2
    continuous = todo
    
    
    if (aggregationFactor > 1) {
      
      print(paste("aggregating by factor of",aggregationFactor))
      continuous <- aggregate(continuous, fact=aggregationFactor, fun=mean)
    } else {
      print("skipping aggregation")
      
    }
    
    plot(continuous[[1]])
    points(sample.locales)
    dev.off()
    
    ## need to transform the abundances by Inverse Monomolecular (7)
    r.tran <- Resistance.tran(transformation = 7,
                              shape = 3.5,
                              max = 150,
                              r = continuous)
    
    values(r.tran) = rescale(values(r.tran))
    
    png(paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".FULL_RESISTANCE_SURFACE.png",sep=""))
    plot(r.tran)
    dev.off()
    
    writeRaster(r.tran,paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".FULL_RESISTANCE_SURFACE.asc",sep=""),
                format="ascii")
    
    
    
    #print("making temps")
    temp = sample.locales
    temp = SpatialPoints(temp)
    #temp = temp[!is.na(extract(continuous,temp)),]
    temp = temp[complete.cases(extract(continuous,temp)),]
    
    individuals = as.numeric(rownames(temp@coords))
    included = samples[individuals,]
    #catalog = included$CATALOG.NUMBER
    latlong = paste(included$lat,included$long)
    
    ## get euclidian distances
    
    rawvalues = as.data.frame(extract(continuous,temp))
    rownames(rawvalues) = latlong
    rawdist = as.data.frame(as.matrix(dist(rawvalues)))
    rawfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".RAW-EUCLIDIAN.txt",sep="")
    write.csv(rawdist,rawfile)
    
    
    outfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".MORPH-AND-GENE-DISTANCES.txt",sep="")
    #print(paste("outfile:",outfile))
    #print(outfile)
    
    #print("gdist prep")
    gdist.inputs <- gdist.prep(length(temp),
                               samples = temp,
                               ##longlat = F, ## may need to comment
                               method = 'costDistance') ## works with costdist?
    
    print("gdist response")
    
    gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                    #r = continuous,
                                    #r=Resist,
                                    r=r.tran,
                                    scl=F) # works if latlongF and methodCost
    
    #print("matrix")
    distdf = as.matrix(gdist.response)
    #print("rename")
    rownames(distdf) = latlong
    colnames(distdf) = latlong
    
    #print("output")
    write.csv(distdf,outfile)
    new <- Sys.time() - old # calculate difference
    print(paste("Time elapsed:",new))
    
  }
}

if (dosmallnospp == T) {
  
  samples = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/LATLONGS_ONLY.txt",sep="\t")
  
  whereraster="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/"
  #whereraster="/home/kprovost/nas1/ENMs/RESCALED/"
  
  pathlist = list.files(path=whereraster,
                        pattern=glob2rx(paste("*rescaled*.asc",sep="")),
                        recursive = T)
  pathlist = pathlist[lapply(pathlist,function(x) length(grep("USA",x,value=FALSE))) == 0]
  pathlist = pathlist[lapply(pathlist,function(x) length(grep("zoomed",x,value=FALSE))) == 0]
  
  sample.locales = (samples[,c(2,1)])
  sample.locales = unique(sample.locales)
  sample.locales = sample.locales[complete.cases(sample.locales),]
  
  maxlat = max(sample.locales$lat,na.rm = T)+1
  maxlong = max(sample.locales$long,na.rm = T)+1
  minlat = min(sample.locales$lat,na.rm = T)-1
  minlong = min(sample.locales$long,na.rm = T)-1
  
  print(paste(minlong,maxlong,minlat,maxlat))
  
  print("setting extent")
  ext = extent(minlong,maxlong,minlat,maxlat)
  print(ext)
  
  sample.locales = SpatialPoints(sample.locales)
  
  for (i in 1:length(pathlist)) {
    
    print(paste(i,"of",length(pathlist)))
    old <- Sys.time()
    print(old)
    
    rasterpath = paste(whereraster,pathlist[i],sep="")
    
    print(rasterpath)
    
    splitpath = strsplit(rasterpath,"/")[[1]]
    suffix = splitpath[length(splitpath)]
    
    #print("load raster")
    continuous = raster(rasterpath)
    
    ## crop raster
    print("cropping")
    continuous = crop(continuous,ext)
    png(paste(outputdirectory,suffix,".AGGREGATE-",aggregationFactor,".png",sep=""))
    par(mfrow=c(1,2))
    plot(continuous)
    
    
    if (aggregationFactor > 1) {
      
      print(paste("aggregating by factor of",aggregationFactor))
      continuous <- aggregate(continuous, fact=aggregationFactor, fun=mean)
    } else {
      print("skipping aggregation")
      
    }
    
    plot(continuous)
    points(sample.locales)
    dev.off()
    
    temp = sample.locales
    temp = SpatialPoints(temp)
    temp = temp[complete.cases(extract(continuous,temp)),]
    
    individuals = as.numeric(rownames(temp@coords))
    included = samples[individuals,]
    latlong = paste(included$lat,included$long)
    
    rawvalues = as.data.frame(extract(continuous,temp))
    rownames(rawvalues) = latlong
    rawdist = as.data.frame(as.matrix(dist(rawvalues)))
    rawfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".RAW-EUCLIDIAN.txt",sep="")
    write.csv(rawdist,rawfile)
    
    
    outfile = paste(outputdirectory,suffix,"AGGFACT-",aggregationFactor,".MORPH-AND-GENE-DISTANCES.txt",sep="")
    
    gdist.inputs <- gdist.prep(length(temp),
                               samples = temp,
                               method = 'costDistance') ## works with costdist?
    
    print("gdist response")
    gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
                                    r = continuous,
                                    scl=F) # works if latlongF and methodCost
    
    distdf = as.matrix(gdist.response)
    rownames(distdf) = latlong
    colnames(distdf) = latlong
    
    write.csv(distdf,outfile)
    new <- Sys.time() - old # calculate difference
    print(paste("Time elapsed:",new))
    
  }
}



