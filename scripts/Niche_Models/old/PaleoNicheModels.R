require(dismo)
require(ENMeval)
require(phyloclim)
require(sp)
require(rgdal)
require(rgeos)
require(ENMTools)
require(spocc)
taxon = 'Ambystoma_tigrinum'

#Build a best model
#get occ
occ = occ(taxon, from = 'gbif', limit =900)
occdf = occ2df(occ)
summary(occ)

#get climate
# Env = stack("/data/spbio/climgrids/bio.gri")
Env = stack("~/envirem.gri")
ext = extent(c(-125, -40, 25, 60))
Env = crop(Env, ext)

#clean
loc = cbind(occdf$longitude, occdf$latitude)
loc = loc[loc[,1]<= -40,]
loc = na.omit(loc)
df = data.frame(occdf)
extr = extract(Env, cbind(df$longitude, df$latitude))
df = df[!is.na(extr[,1]),]

#thin
require(spThin)
df = data.frame(occdf)
extr = extract(Env, cbind(df$longitude, df$latitude))
df = df[!is.na(extr[,1]),]

thin<-thin(loc.data = df, 
               lat.col = "latitude", 
               long.col = "longitude",
               spec.col = "name", 
               thin.par = 10, 
               reps = 10, 
               locs.thinned.list.return = T, 
               write.files = T, 
               max.files = 2, 
               out.dir = paste(taxon, "/", sep=''), 
               write.log.file = T)

newthin = read.csv(paste(taxon, "/thinned_data_thin1.csv", sep=''))

#Buffer background sampling: 200km
#library(ENMTools) #if not already loaded
#with all points (before thinning)
#bg1 = background.points.buffer(df[,2:3], radius = 500000, n = 5000, mask = Env[[1]])
#or after thinning
bg2 = background.points.buffer(newthin[,2:3], radius = 200000, n = 5000, mask = Env[[1]])

plot(Env[[1]], col = viridis::viridis(99))
points(bg2, col = 'grey')
points(newthin[,2:3], col = 'red', pch = 20)


#ENMeval model testing
thinres = ENMevaluate(occ=newthin[,2:3], env = Env, method='block', parallel=T, numCores=4, fc=c("L", "LQ", "H"), RMvalues=seq(1,4,1), rasterPreds=F)
print(thinres@results)
eval.plot(thinres@results, "Mean.AUC", )


#predictbest
setsort = thinres@results[order(thinres@results[,'Mean.ORmin']),]
setsort2 = setsort[order(setsort[,'Mean.AUC'], decreasing=TRUE),]
top = setsort2[1,]
best.thin = which(as.character(thinres@results[,1]) == as.character(setsort2[1,1]))
pred.thin = predict(Env, thinres@models[[best.thin]])
plot(pred.thin, col=viridis::viridis(99))

#thresholding
ev.set <- evaluate(newthin[,2:3], thinres@bg.pts, thinres@models[[best.thin]], Env)
th1 = threshold(ev.set)
p1.nomit = pred.thin >= th1$no_omission
p1.equal = pred.thin >= th1$equal_sens_spec ##probably this one is better
Paleoclimate
lgm = stack('/data/spbio/climgrids/miroc/lgm_envirem.gri')
names(lgm) = names(Env) #Make names match 
lgm = crop(lgm, extent(Env))
lgm.pred.thin = predict(lgm,  thinres@models[[best.thin]]) 

lgm.binary = lgm.pred.thin>= th1$equal_sens_spec

plot(lgm.binary+p1.equal, col = c('black','red', 'green', 'blue'))