## niche comparisons 
library(spocc)
library(raster)
library(viridis)
library(dismo)
library(ecospat)
# devtools::install_github("cran/ecospat")
library(ENMTools)

## first download data for species
d1 = occ('Toxostoma crissale')
d2 = occ('Toxostoma lecontei')
occd1 = data.frame(occ2df(d1))
occd2 = data.frame(occ2df(d2))

## get worldclim
bio = getData("worldclim",var="bio",res=2.5) ## smaller is 0.5
## worthwhile to set a directory where you're keeping the data
## rather than downloading all the time
## you can also download=F and then save them manually as TIF or ASC
sta = bio

## crop worldclim
## set extent to not go into atlantic ocean
occd1 = occd1[occd1$lon <= -40,]
occd2 = occd2[occd2$lon <= -40,]
occd1 = na.omit(occd1)
occd2 = na.omit(occd2)
## then clip to data
combine.lat = c(occd1$lat, occd2$lat)
combine.lon = c(occd1$lon, occd2$lon)
ext=extent(c(min(combine.lon)-5, max(combine.lon)+5, min(combine.lat)-5, max(combine.lat)+5))
Env = crop(sta, ext)

## see the data 
png("Toxostoma_crissale-lecontei_occurrences.png")
plot(Env[[1]], col =viridis(99))
points(occd1[,2:3], pch=21, col = 'darkred',bg=rgb(0,0,0,0))
points(occd2[,2:3], pch=22, col = 'grey',bg=rgb(0,0,0,0))
dev.off()

## now use ecospat to do niche similarity
## first dobackground points, buffer around the points
## radius is in m, so this is 200 km
## n is the number of bg points to select, so in this case 100xocc records
bg1 = background.points.buffer(occd1[,2:3], radius = 200000,
                               n = 100*nrow(occd1), mask = Env[[1]])
bg2 = background.points.buffer(occd2[,2:3], radius = 200000, 
                               n = 100*nrow(occd2), mask = Env[[1]])
## pull out env data for those bg points
extract1 = na.omit(cbind(occd1[,2:3], 
                         extract(Env, occd1[,2:3]), rep(1, nrow(occd1))))
extract2 = na.omit(cbind(occd2[,2:3], 
                         extract(Env, occd2[,2:3]), rep(1, nrow(occd2))))
colnames(extract1)[ncol(extract1)] = 'occ'
colnames(extract2)[ncol(extract2)] = 'occ'
## combine the points with the data 
extbg1 = na.omit(cbind(bg1, extract(Env, bg1), rep(0, nrow(bg1))))
extbg2 = na.omit(cbind(bg2, extract(Env, bg2), rep(0, nrow(bg2))))
colnames(extbg1)[ncol(extbg1)] = 'occ'
colnames(extbg2)[ncol(extbg2)] = 'occ'
dat1 = rbind(extract1, extbg1)
dat2 = rbind(extract2, extbg2)

## do a pca of the background points
pca.env <- dudi.pca(
  rbind(dat1, dat2)[,3:21],
  scannf=FALSE,
  nf=2
)
## see the contribution of each variable from the pca
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

## extract the pca values for both species
scores.globclim<-pca.env$li # PCA scores for the whole study area (all points)
scores.sp1 <- suprow(pca.env,
                     extract1[which(extract1[,22]==1),3:21])$li # PCA scores for the species 1 distribution
scores.sp2 <- suprow(pca.env,
                     extract2[which(extract2[,22]==1),3:21])$li # PCA scores for the species 2 distribution
scores.clim1 <- suprow(pca.env,dat1[,3:21])$li # PCA scores for the whole native study area species 1
scores.clim2 <- suprow(pca.env,dat2[,3:21])$li # PCA scores for the whole native study area species 2

## make a dynamic occurrence densities grid
## grid of occ densities along one or two environmental gradients
## glob = env variables in background 
## glob1 = env variables for species
## sp = occurrences of species
## R = resolution
## th.sp = a threshhold to elimite low density values of species occurrences
grid.clim1 <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim1,
  sp = scores.sp1,
  R = 100,
  th.sp = 0
)
grid.clim2 <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim2,
  sp = scores.sp2,
  R = 100,
  th.sp = 0
)

## is this a projection of the niche?
## its the pca space occupied by the two species
png("Toxostoma_pcaSpace_cris-leco.png")
par(mfrow=c(1,2))
plot(grid.clim1$w,main="Toxostoma crissale")
plot(grid.clim2$w,main="Toxostoma lecontei")
dev.off()

## calculate overlap between the species niches
## Schoener's D
D.overlap <- ecospat.niche.overlap (grid.clim1, grid.clim2, cor=T)$D 
D.overlap # 0.3298541 for Toxostoma cris vs leco
## modified Hellinger's I
I.overlap = ecospat.niche.overlap (grid.clim1, grid.clim2, cor=T)$I
I.overlap # 0.5620279 for Toxostoma cris vs leco
## this is looking for overlap in PCA-environmental space of the two species

## test for niche equivalence
## this is the observed niche overlap when you randomize species id
## greater means you test for niche conservatism
## lower means you test for niche divergence 
eq.test <- ecospat.niche.equivalency.test(grid.clim1, grid.clim2,
                                          rep=10, alternative = "greater") ##rep = 1000 recommended for operational runs
## for 10 takes 1 minute

## then test for niche similarity -- 
## overlap between spp 1 and overlaps between random spp 2 bg niches
## rand.type 2 means only z2 is randomly shifted
sim.test <- ecospat.niche.similarity.test(grid.clim1, grid.clim2,
                                          rep=1000, alternative = "greater",
                                          rand.type=2) 
sim.test2 <- ecospat.niche.similarity.test(grid.clim2,grid.clim1,
                                          rep=1000, alternative = "greater",
                                          rand.type=2) 
png("Toxostoma_D-equivalency_cris-leco.png")
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()
png("Toxostoma_D-similarity_cris-leco.png")
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()
png("Toxostoma_D-similarity_leco-cris.png")
ecospat.plot.overlap.test(sim.test2, "D", "Similarity")
dev.off()

## try running with lesser
eq.test1 <- ecospat.niche.equivalency.test(grid.clim1, grid.clim2,
                                          rep=10, alternative = "lesser") ##rep = 1000 recommended for operational runs
sim.test3 <- ecospat.niche.similarity.test(grid.clim1, grid.clim2,
                                          rep=1000, alternative = "lesser",
                                          rand.type=2) 
sim.test4 <- ecospat.niche.similarity.test(grid.clim2,grid.clim1,
                                           rep=1000, alternative = "lesser",
                                           rand.type=2) 
png("Toxostoma_D-equivalency-lesser_cris-leco.png")
ecospat.plot.overlap.test(eq.test1, "D", "Equivalency")
dev.off()
png("Toxostoma_D-similarity-lesser_cris-leco.png")
ecospat.plot.overlap.test(sim.test3, "D", "Similarity")
dev.off()
png("Toxostoma_D-similarity-lesser_leco-cris.png")
ecospat.plot.overlap.test(sim.test4, "D", "Similarity")
dev.off()