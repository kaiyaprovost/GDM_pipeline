#detach("package:subsppLabelR", unload = TRUE)
#library(devtools)
#devtools::install_github('kaiyaprovost/subsppLabelR', force = F)
# devtools::install_github("cran/ecospat")
# devtools::install_github("danlwarren/ENMTools")
# devtools::install_github("rsh249/rasterExtras")
# devtools::install_github("petrelharp/landsim")
library(subsppLabelR)
library(ecospat)
library(ENMTools)
library(rasterExtras)
library(landsim)
library(sys)
library(scales)

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

path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/"
setwd(path)

specieslist=c(
  "Amphispiza bilineata",
  "Auriparus flaviceps",
  "Campylorhynchus brunneicapillus",
  "Cardinalis sinuatus",
  #"Melozone crissalis",
  "Melozone fusca",
  #"Melozone fusca nominate",
  "Phainopepla nitens",
  "Polioptila melanura",
  "Toxostoma crissale",
  "Toxostoma curvirostre",
  #"Toxostoma curvirostre nominate",
  #"Toxostoma palmeri",
  "Vireo bellii"
)

## stack time slices
for (spp in rev(specieslist)) {
  
  print(spp)
  current_f = list.files(path=path,pattern=paste("Thresh_EqualSensSpec_",spp,"_worldclim.asc$",sep=".+"),recursive = T,full.names = T)[1]
  current = raster (current_f)
  
  LGM_f = list.files(path=path,pattern=paste("Thresh_EqualSensSpec_",spp,"_LGM.asc$",sep=".+"),recursive = T,full.names = T)[1]
  LGM = raster (LGM_f)
  
  MID_f = list.files(path=path,pattern=paste("Thresh_EqualSensSpec_",spp,"_MID.asc$",sep=".+"),recursive = T,full.names = T)[1]
  MID = raster (MID_f)
  
  #current = raster(paste("Thresh_EqualSensSpec_BestModel_",spp,"_addedLayers_zoomedin.asc",sep=""))
  #LGM = raster(paste("Thresh_EqualSensSpec_BestModel_",spp,"_addedLayers_LGM.asc",sep=""))
  #MID = raster(paste("Thresh_EqualSensSpec_BestModel_",spp,"_addedLayers_MID.asc",sep=""))
  # stacked = stack(list(current,LGM,MID))
  
  
  
  # png(paste("Thresh_Timestack_",spp,".png",sep=""))
  # plotRGB(stacked, r=1, g=2, b=3,scale=1,colNA="grey")
  # legend(legend=c("Now-Only","LGM-Only","Mid-Only",
  #                 "Now-LGM","Now-Mid","LGM-Mid",
  #                 "None","All"),
  #        x="right",
  #        title=spp,
  #        bty="n",
  #         fill=c(rgb(1,0,0),rgb(0,1,0),rgb(0,0,1),
  #                rgb(1,1,0),rgb(1,0,1),rgb(0,1,1),
  #                rgb(0,0,0),rgb(1,1,1))
  #        #fill=c("#E31A1C","#33A02C","#1F78B4",
  #       #        "#FB9A99","#6A3D9A","#A6CEE3",
  #       #        "black","white")
  #        )
  # dev.off()
  
  MID_2 = MID
  values(MID_2)[values(MID_2)==1 & !(is.na(values(MID_2)))] = 2
  LGM_4 = LGM
  values(LGM_4)[values(LGM_4)==1 & !(is.na(values(LGM_4)))] = 4
  
  summed = current+MID_2+LGM_4
  
  png(paste("Thresh_Timestack_",spp,"_difcolors.png",sep=""))
  plot(summed,col=c(
    "black", #//
    "#FF7F00", #C
    "#FB9A99", #M
    "#E31A1C", #CM
    "#1F78B4", #L
    "#6A3D9A", #CL
    "#CAB2D6", #ML
    "grey" #CLM
    ))
  legend(legend=c("None","Now-Only","Mid-Only","Now-Mid",
                  "LGM-Only","Now-LGM","Mid-LGM","All"),
         x="right",
         title=spp,
         bty="n",
         fill=c("black", #//
                "#FF7F00", #C
                "#FB9A99", #M
                "#E31A1C", #CM
                "#1F78B4", #L
                "#6A3D9A", #CL
                "#CAB2D6", #ML
           "grey" #CLM
           ))
  dev.off()
  
  
  # current2 = raster(paste("AddedPredThin_BestModel_",spp,"_addedLayers_zoomedin.asc",sep=""))
  # LGM2 = raster(paste("AddedPredThin_BestModel_",spp,"_addedLayers_LGM.asc",sep=""))
  # MID2 = raster(paste("AddedPredThin_BestModel_",spp,"_addedLayers_MID.asc",sep=""))
  # values(current2) = scales::rescale(values(current2),to=c(0,1))
  # values(LGM2) = scales::rescale(values(LGM2),to=c(0,1))
  # values(MID2) = scales::rescale(values(MID2),to=c(0,1))
  # stacked2 = stack(list(current2,LGM2,MID2))
  # png(paste("Continuous_Timestack_",spp,".png",sep=""))
  # plotRGB(stacked2, r=1, g=2, b=3,scale=1,colNA="grey",stretch="hist",
  #         xlab=spp)
  # legend(legend=c("Now-Only","LGM-Only","Mid-Only",
  #                 "Now-LGM","Now-Mid","LGM-Mid",
  #                 "None","All"),
  #        title=spp,
  #        x="right",
  #        bty="n",
  #        # fill=c(rgb(1,0,0),rgb(0,1,0),rgb(0,0,1),
  #        #        rgb(1,1,0),rgb(1,0,1),rgb(0,1,1),
  #        #        rgb(0,0,0),rgb(1,1,1))
  #        fill=c("#E31A1C","#33A02C","#1F78B4",
  #               "#FB9A99","#6A3D9A","#A6CEE3",
  #               "black","white")
  #        )
  # dev.off()
  
}


path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WITH_MEXICO/"
setwd(path)
files = list.files(pattern="rescaled.asc$",recursive=T)

numbins = 10
binsize = 1/numbins

for (asc in files) {
  print(asc)
  print("--reading")
  ras = raster(asc)
  spp = strsplit(strsplit(asc,"/")[[1]][3],"_")[[1]][3]
  floorras = ras
  png(paste(asc,"_",numbins,"bins.png",sep=""))
  print("--converting")
  values(floorras) = round(values(ras)/binsize)*binsize
  plot(floorras, col=viridis::viridis(99),main=spp)
  dev.off()
  
}


## see how USA vs worldclim rasters compare

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/")
current = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/USA_ONLY/Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_USA.asc_-5.90643692016602_-13.0102529525757.rescaled.asc"
past = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc_-6.46640157699585_-69.5545196533203.rescaled.asc"

curras = raster(current)
pasras = raster(past)
pasras = crop(pasras,extent(curras))

extent(curras)
extent(pasras)

#paslrg = resample(pasras,curras)
#extent(paslrg)

cursml = resample(curras,pasras)
values(cursml) = scales::rescale(values(cursml),to=c(0,1))
values(pasras) = scales::rescale(values(pasras),to=c(0,1))


#difras = curras-paslrg
difras = cursml-pasras
difras = crop(difras,extent(curras))
plot(difras)

dif = stack(cursml,pasras)
names(dif) = c("cur","pas")
vals = values(dif)
vals = unique(vals)
vals = vals[complete.cases(vals),]
nrow(vals)

hist(values(difras))

png("test.png")
plot(vals,col=rgb(0,0,0,0.1),pch=16)
abline(0,1,col="red",lwd=2,lty=2)
mod = lm(vals[,2]~vals[,1])
#summary(mod)
abline(mod,col="cyan",lwd=2)
text(x=0.8,y=0.4,labels=c(paste("R^2:",round(summary(mod)$adj.r.squared,2))))
dev.off()

resid = residuals(mod)
hist(resid)
