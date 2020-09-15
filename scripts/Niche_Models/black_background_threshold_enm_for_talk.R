library(raster)

m_p = "/Users/kprovost/Pictures/Thresh_EqualSensSpec_BestModel_Melozone fusca_addedLayers_worldclim.asc"
m_l = "/Users/kprovost/Pictures/Thresh_EqualSensSpec_BestModel_Melozone fusca_addedLayers_LGM.asc"
a_p = "/Users/kprovost/Pictures/Thresh_MaxSpecSens_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc" ## OH NO
a_l = "/Users/kprovost/Pictures/Thresh_MaxSpecSens_BestModel_Amphispiza bilineata_addedLayers_LGM.asc"
a_p_f = "/Users/kprovost/Pictures/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc"
a_l_f = "/Users/kprovost/Pictures/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_LGM.asc"

m_p = raster(m_p)
m_l = raster(m_l)
a_p = raster(a_p)
a_l = raster(a_l)
a_p_f = raster(a_p_f)
a_l_f = raster(a_l_f)

m_p= raster::crop(m_p,extent(c(-124,-93,20,40)))
m_l= raster::crop(m_l,extent(c(-124,-93,20,40)))
a_p= raster::crop(a_p,extent(c(-124,-93,20,40)))
a_l= raster::crop(a_l,extent(c(-124,-93,20,40)))
a_p_f= raster::crop(a_p_f,extent(c(-124,-93,20,40)))
a_l_f= raster::crop(a_l_f,extent(c(-124,-93,20,40)))

a_l_f=a_l_f-min(values(a_p_f),na.rm=T)
a_p_f=a_p_f-min(values(a_p_f),na.rm=T)

path <- rgdal::readOGR(dsn = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/WCD_from_Ed_1236_4iouz24kv_ConDivides.gpx", layer = "tracks")
path= raster::crop(path,extent(c(-124,-93,20,40)))

pdf("bil_magenta_suitability.pdf",height=2.5,width=6)
par(mfrow=c(1,2),mar=c(0,0,0,0))
x=colorRampPalette(c("hotpink4","#E7298A","pink","moccasin","yellow","white"))
plot((a_p_f),col=rev(x(10)),zlim=c(14,32.3),colNA=rgb(0,0,0,0),legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
plot((a_l_f),col=rev(x(10)),zlim=c(14,32.3),colNA=rgb(0,0,0,0),legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
dev.off()

pdf("bil_magenta_thresh.pdf",height=2.5,width=6)
par(mfrow=c(1,2),mar=c(0,0,0,0))
x=colorRampPalette(c("hotpink4","#E7298A","pink","yellow","white"))
plot((a_p),col=c("grey","#E7298A"),colNA=rgb(0,0,0,0),legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
plot((a_l),col=c("grey","#E7298A"),colNA=rgb(0,0,0,0),legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
dev.off()

m_p[values(m_p)==1] = 2
a_p[values(a_p)==1] = 2
m = m_p+m_l
a = a_p+a_l


orange="#E7298A"
red="violetred4"
yellow="magenta4"

pdf("bil_fus_lgm_present_timestack_only.pdf",height=2.5,width=6)
par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(a,ylim=c(20,40),xlim=c(-124,-93),col=c("grey",yellow,red,orange),
     colNA="black",legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
plot(m,ylim=c(20,40),xlim=c(-124,-93),col=c("grey",yellow,red,orange),
     colNA="black",legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
dev.off()

## timestacks for all ten species 
lgms = list.files(path="/Users/kprovost/Downloads/climate thresh/",
                  pattern="LGM.asc",full.names = T)
pres = list.files(path="/Users/kprovost/Downloads/climate thresh/",
                  pattern="worldclim.asc",full.names = T)
path <- rgdal::readOGR(dsn = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/WCD_from_Ed_1236_4iouz24kv_ConDivides.gpx", layer = "tracks")
path= raster::crop(path,extent(c(-124,-93,20,40)))
setwd("~")
for(i in 1:length(lgms)){
  p = raster(pres[i])
  l = raster(lgms[i])
  p = raster::crop(p,extent(c(-124,-93,20,40)))
  l = raster::crop(l,extent(c(-124,-93,20,40)))
  
  pdf(paste(basename(pres[i]),"_magenta_thresh.pdf"),height=2.5,width=3)
  par(mar=c(0,0,0,0))
  x=colorRampPalette(c("hotpink4","#E7298A","pink","yellow","white"))
  plot((p),col=c("grey","#E7298A"),colNA=rgb(0,0,0,0),legend=F,yaxt="n",xaxt="n")
  plot(path,add=T,col="black")
  dev.off()
  pdf(paste(basename(lgms[i]),"_magenta_thresh.pdf"),height=2.5,width=3)
  par(mar=c(0,0,0,0))
  plot((l),col=c("grey","#E7298A"),colNA=rgb(0,0,0,0),legend=F,yaxt="n",xaxt="n")
  plot(path,add=T,col="black")
  dev.off()
  
  p[values(p)==1] = 2
  m = p+l
  orange="#E7298A"
  red="violetred4"
  yellow="magenta4"
  
  pdf(paste(basename(lgms[i]),"_magenta_STACK.pdf"),height=2.5,width=3)
  par(mar=c(0,0,0,0))
  x=colorRampPalette(c("hotpink4","#E7298A","pink","yellow","white"))
  plot(m,ylim=c(20,40),xlim=c(-124,-93),col=c("grey",yellow,red,orange),
       colNA="black",legend=F,yaxt="n",xaxt="n")
  plot(path,add=T,col="black")
  dev.off()
  
}

## abundance

abun_a = raster("/Users/kprovost/Downloads/clipped/BILINEATA_Idw_abundance_clipped.asc")
abun_m = raster("/Users/kprovost/Downloads/clipped/FUSCA_Idw_abundance_clipped.asc")

breaks=c(0,1,3,10,30,100)
x=colorRampPalette(c("#1B9E77","white"))
x(10)

pdf("bil_fus_abun.pdf",height=2.5,width=6)
par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(abun_a,ylim=c(20,40),xlim=c(-124,-93),#breaks=breaks,
     col=rev(x(7)[1:6]),
     colNA="black",legend=T,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
plot(abun_m,ylim=c(20,40),xlim=c(-124,-93),breaks=breaks,
     col=rev(x(7)[1:6]),
     colNA="black",legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
dev.off()

## abun all species
path <- rgdal::readOGR(dsn = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/WCD_from_Ed_1236_4iouz24kv_ConDivides.gpx", layer = "tracks")
path= raster::crop(path,extent(c(-124,-93,20,40)))
mex = shapefile("/Users/kprovost/Downloads/enm_layers/mexstates/mexstates.shp")
usa = shapefile("/Users/kprovost/Downloads/enm_layers/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
mex = raster::crop(mex,extent(c(-124,-93,20,40)))
usa=raster::crop(usa,extent(c(-124,-93,20,40)))

setwd("~")
for (file in list.files(path="/Users/kprovost/Downloads/clipped abun/",full.names = T)){
  breaks=c(0,1,3,10,30,100)
  x=colorRampPalette(c("#1B9E77","white"))
  #x=colorRampPalette(c("#d7191c","#fdae61","#ffffbf",
  #  "#abdda4","#2b83ba"))
  #x(10)
  abun_a = raster(file)
  abun_a2 <- crop(abun_a, extent(usa))
  abun_a3 <- mask(abun_a2, usa)
  
  pdf(paste(basename(file),"log 3.pdf"),height=2.5,width=3)
  par(mar=c(0,0,0,0))
  plot((abun_a3),
       #ylim=c(20,40),xlim=c(-124,-93),#breaks=breaks,
       ylim=c(25,40),xlim=c(-124,-93),#breaks=breaks,
       col="grey",zlim=c(0,100),
       colNA="black",legend=F,yaxt="n",xaxt="n")
  plot(log10(abun_a3),
       #ylim=c(20,40),xlim=c(-124,-93),#breaks=breaks,
       ylim=c(25,40),xlim=c(-124,-93),#breaks=breaks,
       col=rev(x(7)[1:6]),zlim=c(0,2),add=T,
       #col=rev(x(100)),zlim=c(0,2),add=T,
       legend=T,yaxt="n",xaxt="n")
  #plot(path,add=T,col="black")
  #plot(mex,add=T,density=25,col="black")
  #plot(usa,add=T,density=0,col="black")
  dev.off()
}


## bioclim
brown="#A6761D"
yellow="#E6AB02"
x=colorRampPalette(c(brown,"lightyellow"))

temp=raster("/Users/kprovost/Downloads/bio_2-5m_bil/bio1.bil")
temp=raster::crop(temp,extent(c(-124,-93,20,40)))
prec=raster("/Users/kprovost/Downloads/bio_2-5m_bil/bio12.bil")
prec=raster::crop(prec,extent(c(-124,-93,20,40)))

pdf("bioclim_temp_logprec_brown.pdf",height=2.5,width=6)
par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(temp,ylim=c(20,40),xlim=c(-124,-93),
     col=x(10),
     colNA="black",legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
plot(log(prec),ylim=c(20,40),xlim=c(-124,-93),
     col=rev(x(10)),
     colNA="black",legend=F,yaxt="n",xaxt="n")
plot(path,add=T,col="black")
dev.off()

x="#7570b3ff"
y="#e6c0ffff"
