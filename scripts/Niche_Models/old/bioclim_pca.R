## pcas of bioclim variables plus elevation

library(raster)

folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/bio_2-5m_bil/"
worldclim = list.files(folder,pattern="bil$",full.names = T)
clim = raster::stack(worldclim)
## clip to SW
#ext = raster::extent(c(-125,-99, 22, 44))
#ext = raster::extent(c(-118,-98, 21, 36))
ext = raster::extent(c(-125,-69,10,55))
elev=raster::raster("/Users/kprovost/Documents/ENMs/OLD/NASA_TOPOGRAPHY_SRTM_RAMP2_TOPO_2000-02-11_gs_3600x1800.asc")

clim = raster::crop(clim, ext)
elev = raster::crop(elev, ext)

plot(clim[[1]])
plot(elev)

bigclim = raster::raster("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer4.tif")


# longs=c(-101.00,-109.03965)
# lats = c(25.67,32.0526078)
# x = coordinates((cbind(longs,lats)))
# y = extract(clim[[c("bio1","bio12")]],x)
# png("precip_barplot.png",height=280)
# par(mar=c(2,4,0,0))
# barplot(y[,2],horiz=F,las=1,col=c("lightblue","blue"),ylim=c(0,310),names=c("Coahuila","New Mexico"),ylab="Precipitation (mm)")
# dev.off()
# png("temp_barplot.png",height=280)
# par(mar=c(2,4,0,0))
# barplot(y[,1]/10,horiz=F,las=1,col=c("red","pink"),names=c("Coahuila","New Mexico"),ylim=c(0,21),ylab="Temperature (°C)")
# dev.off()

shapefile = raster::shapefile("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/deserts_north_america_wwf_SONCHI.shp")
plot(shapefile,add=F)
clim = raster::crop(clim,shapefile)
plot(clim[[1]])
plot(shapefile,add=T)
elev = raster::crop(elev,shapefile)
plot(elev)
plot(shapefile,add=T)

clim2 = clim
r.polys <- rasterize(shapefile, clim2, field = 1, fun = "mean", 
                     update = F, updateValue = "NA")
plot(r.polys)
clim2[is.na(r.polys),] = NA
plot(clim2[[1]])


plot(clim2[[c("bio1","bio4","bio12","bio15")]])



elev2 = elev
r.polys <- rasterize(shapefile, elev2, field = 1, fun = "mean", 
                     update = F, updateValue = "NA")
plot(r.polys)
elev2[is.na(r.polys),] = NA
plot(elev2)

data = values(clim2)
data2 = unique(data)
data2 = data2[complete.cases(data2),]
data3 = data2
data3[,c(1:3,12:19)] = data3[,c(1:3,12:19)]/10

dataelev = values(elev2)
dataelev2 = unique(dataelev)
dataelev2 = dataelev2[complete.cases(dataelev2)]

climr = resample(clim2[[1]],elev2)
climp = resample(clim2[[4]],elev2)

## correlate values
dat1 = values(climr)
dat2 = values(climp)
dat = cbind(dataelev,dat1,dat2)
dat=dat[complete.cases(dat),]
par(mfrow=c(1,2))
plot(dat[,1],dat[,2])
abline(lm(dat[,2]~dat[,1]),col="red")
plot(dat[,1],dat[,3])
abline(lm(dat[,3]~dat[,1]),col="blue")


sds=matrixStats::colSds(data3)
means=colMeans(data3)
mins = sapply(1:ncol(data3),FUN = function(x){min(data3[,x])})
maxs = sapply(1:ncol(data3),FUN = function(x){max(data3[,x])})

summaries=rbind(sds,means,mins,maxs)
summaries

pca = prcomp(data2,center = T,scale. = T,rank. = 3)
#newdata = pca$x
pca
abs(pca$rotation)
corrplot::corrplot(abs(t(pca$rotation)),is.corr=F,cl.lim=c(min(abs(pca$rotation)),
                                                           max(abs(pca$rotation))),
                   method="number")
summary(pca)
# r=3,g=1,b=2

corrplot::corrplot(abs(t(pca$rotation)),is.corr=F,cl.lim=c(min(abs(pca$rotation)),
                                                           max(abs(pca$rotation))),
                   method="ellipse",
                   col=c(rev(viridis::viridis(50)),viridis::viridis(50)))


pred <- predict(pca, newdata=values(clim))
pred2 = pred

pred2[,1] = scales::rescale(pred[,1],to=c(0,1))
pred2[,2] = scales::rescale(pred[,2],to=c(0,1))
pred2[,3] = scales::rescale(pred[,3],to=c(0,1))

pred3 <- predict(pca, newdata=data)
pred4 = pred3
pred4[,1] = scales::rescale(pred3[,1],to=c(0,1))
pred4[,2] = scales::rescale(pred3[,2],to=c(0,1))
pred4[,3] = scales::rescale(pred3[,3],to=c(0,1))




climpca = clim[[1:3]]
values(climpca) = pred2
names(climpca)

writeRaster(climpca,"Deserts_Bioclim_PCA_SONCHI.tif",format="GTiff")

#png("Deserts_Bioclim_PCA_CONTINENT.png")
png("Deserts_Bioclim_PCA_SONCHI.png")
plotRGB(climpca,scale=1,r=2,g=1,b=3,colNA="white") #312 is best previously
plot(shapefile,add=T)
#points(x,col="black",cex=2,pch=4,lwd=5)
dev.off()

png("Deserts_Bioclim_PCA_SONCHI_ONLY.png")
climpca[is.na(r.polys),] = NA
plotRGB(climpca,scale=1,r=2,g=1,b=3,colNA="white") #312 is best previously
dev.off()


climpca2 = clim[[1:3]]
values(climpca2) = pred4
names(climpca2)
plotRGB(climpca2,scale=1,r=3,g=1,b=2) #1
plot(shapefile,add=T)

plotRGB(climpca,scale=1,r=1,g=1,b=1,colNA="magenta") #1
plotRGB(climpca,scale=1,r=2,g=2,b=2,colNA="magenta") #1
plotRGB(climpca,scale=1,r=3,g=3,b=3,colNA="magenta") #1

#par(mfrow=c(1,3))
#plotRGB(climpca2,scale=1,r=1,g=1,b=1,colNA="magenta") #1
#plotRGB(climpca2,scale=1,r=2,g=2,b=2,colNA="magenta") #1
#plotRGB(climpca2,scale=1,r=3,g=3,b=3,colNA="magenta") #1

plot(climpca2)

temp = clim[[c(1,10+5,10+6,10+7,11-8)]]
names(temp) = toupper(names(temp))
prec = clim[[c(12-8,13-8,14-8,18-8,19-8)]]
names(prec) = toupper(names(prec))

temp_pca = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/temperature_pca.csv",row.names = 1)
temp_pca = temp_pca[1:5,]
prec_pca = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/precipitation_pca.csv",row.names = 1)
prec_pca = prec_pca[1:5,]

temp_val = values(temp)
temp_val = as.data.frame(scale(temp_val, center = TRUE, scale = TRUE))

corrplot::corrplot(cor(temp_val,use="pairwise.complete.obs"),method="number")

temp_pca_val = temp_val
for (i in 1:nrow(temp_pca_val)) { 
  
  if (i %% 100 == 0) {
    print(paste(i/100,"of",nrow(temp_pca_val)/100))
  }
  
  row = t(temp_pca_val[i,])
  newrow = colSums(row * temp_pca)
  temp_pca_val[i,] = newrow
  
} 

temp_pca_ras = temp
values(temp_pca_ras) = as.matrix(temp_pca_val)
writeRaster(temp_pca_ras,"Deserts_Bioclim_TempPCA_From_GDM.asc",format="GTiff")

temp_pca_val[,1] = scales::rescale(temp_pca_val[,1])
temp_pca_val[,2] = scales::rescale(temp_pca_val[,2])
temp_pca_val[,3] = scales::rescale(temp_pca_val[,3])

values(temp_pca_ras) = as.matrix(temp_pca_val)
png("Deserts_Bioclim_TempPCA_From_GDM.png")
plotRGB(temp_pca_ras,scale=1,r=3,g=1,b=2) #1
dev.off()





prec_val = values(prec)
prec_val = as.data.frame(scale(prec_val, center = TRUE, scale = TRUE))

prec_pca_val = prec_val
for (i in 1:nrow(prec_pca_val)) { 
  
  if (i %% 100 == 0) {
    print(paste(i/100,"of",nrow(prec_pca_val)/100))
  }
  
  row = t(prec_pca_val[i,])
  newrow = colSums(row * prec_pca)
  prec_pca_val[i,] = newrow
  
} 

prec_pca_ras = prec
values(prec_pca_ras) = as.matrix(prec_pca_val)
writeRaster(prec_pca_ras,"Deserts_Bioclim_PrecPCA_From_GDM.asc",format="GTiff")

prec_pca_val[,1] = scales::rescale(prec_pca_val[,1])
prec_pca_val[,2] = scales::rescale(prec_pca_val[,2])
prec_pca_val[,3] = scales::rescale(prec_pca_val[,3])

values(prec_pca_ras) = as.matrix(prec_pca_val)
png("Deserts_Bioclim_precPCA_From_GDM.png")
plotRGB(prec_pca_ras,scale=1,r=3,g=1,b=2) #1
dev.off()


prec_pca_val = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/Deserts_Bioclim_PrecPCA_From_GDM.tif")
temp_pca_val = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/Deserts_Bioclim_TempPCA_From_GDM.tif")

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/PrincipalComponents_TempAndPrec.png",
    width=776,height=460)
par(mfrow=c(2,3))
plot(temp_pca_val[[1]],main="PC1T")
plot(temp_pca_val[[2]],main="PC2T")
plot(temp_pca_val[[3]],main="PC3T")
plot(prec_pca_val[[1]],main="PC1P")
plot(prec_pca_val[[2]],main="PC2P")
plot(prec_pca_val[[3]],main="PC3P")
dev.off()




## DO THIS BUT FOR THE TRANSLOCATED ONES
folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/bio_2-5m_bil/"
worldclim = list.files(folder,pattern="bil$",full.names = T)
clim = raster::stack(worldclim)
## clip to SW
ext = raster::extent(c(-125,-99, 22, 44))
#ext = raster::extent(c(-118,-98, 21, 36))

clim = raster::crop(clim, ext)
plot(clim[[1]])

shapefile = raster::shapefile("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/deserts_north_america_wwf.shp")
plot(shapefile,add=T)
clim = raster::crop(clim,shapefile)
plot(clim[[1]])
plot(shapefile,add=T)


clim2 = clim
r.polys <- rasterize(shapefile, clim2, field = 1, fun = "mean", 
                     update = F, updateValue = "NA")
plot(r.polys)

clim2[is.na(r.polys),] = NA

plot(clim2[[1]])

data = values(clim2)
data2 = unique(data)
data2 = data2[complete.cases(data2),]

data3 = data2
data3[,c(1:3,12:19)] = data3[,c(1:3,12:19)]/10

sds=matrixStats::colSds(data3)
means=colMeans(data3)
mins = sapply(1:ncol(data3),FUN = function(x){min(data3[,x])})
maxs = sapply(1:ncol(data3),FUN = function(x){max(data3[,x])})

summaries=rbind(sds,means,mins,maxs)


pca = prcomp(data2,center = T,scale. = T,rank. = 3)
#newdata = pca$x
pca
abs(pca$rotation)
corrplot::corrplot(abs(t(pca$rotation)),is.corr=F,cl.lim=c(min(abs(pca$rotation)),
                                                           max(abs(pca$rotation))),
                   method="number")
summary(pca)
# r=3,g=1,b=2

pred <- predict(pca, newdata=data)
pred2 = pred

pred2[,1] = scales::rescale(pred[,1])
pred2[,2] = scales::rescale(pred[,2])
pred2[,3] = scales::rescale(pred[,3])


climpca = clim[[1:3]]
values(climpca) = pred2
names(climpca)

png("Deserts_Bioclim_PCA.png")
plotRGB(climpca,scale=1,r=3,g=1,b=2) #1
plot(shapefile,add=T)
dev.off()

plotRGB(climpca,scale=1,r=1,g=1,b=1,colNA="magenta") #1
plotRGB(climpca,scale=1,r=2,g=2,b=2,colNA="magenta") #1
plotRGB(climpca,scale=1,r=3,g=3,b=3,colNA="magenta") #1


temp = clim[[c(1,10+5,10+6,10+7,11-8)]]
names(temp) = toupper(names(temp))
prec = clim[[c(12-8,13-8,14-8,18-8,19-8)]]
names(prec) = toupper(names(prec))

temp_pca = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/temperature_pca.csv",row.names = 1)
temp_pca = temp_pca[1:5,]
prec_pca = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/precipitation_pca.csv",row.names = 1)
prec_pca = prec_pca[1:5,]

temp_val = values(temp)
temp_val = as.data.frame(scale(temp_val, center = TRUE, scale = TRUE))

corrplot::corrplot(cor(temp_val,use="pairwise.complete.obs"),method="number")

temp_pca_val = temp_val
for (i in 1:nrow(temp_pca_val)) { 
  
  if (i %% 100 == 0) {
    print(paste(i/100,"of",nrow(temp_pca_val)/100))
  }
  
  row = t(temp_pca_val[i,])
  newrow = colSums(row * temp_pca)
  temp_pca_val[i,] = newrow
  
} 

temp_pca_ras = temp
values(temp_pca_ras) = as.matrix(temp_pca_val)
writeRaster(temp_pca_ras,"Deserts_Bioclim_TempPCA_From_GDM.asc",format="GTiff")

temp_pca_val[,1] = scales::rescale(temp_pca_val[,1])
temp_pca_val[,2] = scales::rescale(temp_pca_val[,2])
temp_pca_val[,3] = scales::rescale(temp_pca_val[,3])

values(temp_pca_ras) = as.matrix(temp_pca_val)
png("Deserts_Bioclim_TempPCA_From_GDM.png")
plotRGB(temp_pca_ras,scale=1,r=3,g=1,b=2) #1
dev.off()





prec_val = values(prec)
prec_val = as.data.frame(scale(prec_val, center = TRUE, scale = TRUE))

prec_pca_val = prec_val
for (i in 1:nrow(prec_pca_val)) { 
  
  if (i %% 100 == 0) {
    print(paste(i/100,"of",nrow(prec_pca_val)/100))
  }
  
  row = t(prec_pca_val[i,])
  newrow = colSums(row * prec_pca)
  prec_pca_val[i,] = newrow
  
} 

prec_pca_ras = prec
values(prec_pca_ras) = as.matrix(prec_pca_val)
writeRaster(prec_pca_ras,"Deserts_Bioclim_PrecPCA_From_GDM.asc",format="GTiff")

prec_pca_val[,1] = scales::rescale(prec_pca_val[,1])
prec_pca_val[,2] = scales::rescale(prec_pca_val[,2])
prec_pca_val[,3] = scales::rescale(prec_pca_val[,3])

values(prec_pca_ras) = as.matrix(prec_pca_val)
png("Deserts_Bioclim_precPCA_From_GDM.png")
plotRGB(prec_pca_ras,scale=1,r=3,g=1,b=2) #1
dev.off()


prec_pca_val = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/Deserts_Bioclim_PrecPCA_From_GDM.tif")
temp_pca_val = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/Deserts_Bioclim_TempPCA_From_GDM.tif")

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/PrincipalComponents_TempAndPrec.png",
    width=776,height=460)
par(mfrow=c(2,3))
plot(temp_pca_val[[1]],main="PC1T")
plot(temp_pca_val[[2]],main="PC2T")
plot(temp_pca_val[[3]],main="PC3T")
plot(prec_pca_val[[1]],main="PC1P")
plot(prec_pca_val[[2]],main="PC2P")
plot(prec_pca_val[[3]],main="PC3P")
dev.off()

