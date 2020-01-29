## RERUN EVERYTHING ONCE DONE TO MAKE SURE GOOD

require(dismo)
require(ENMeval)
require(phyloclim)
require(sp)
require(rgdal)
require(rgeos)
require(ENMTools)
require(spocc)
require(spThin)

## FIRST CARDINALIS SINUATUS
## THEN DID FIND AND REPLACE FOR GLANDULOSA

##Prosopis glandulosa
occ = occ('Prosopis glandulosa')
##occ = occ('Cardinalis sinuatus') ## limit of 500 default
occ
occdf = occ2df(occ) ## this is a ggplot dataframe but its about the same as a normal one
occdf = data.frame(occdf)
summary(occ)

loc = cbind(occdf$longitude, occdf$latitude)
loc = loc[loc[,1]<= -60,] ## this removes everything below -60 degrees longitude
loc = na.omit(loc)
#nrow(loc)
#loc

Env = stack("/data/spbio/climgrids/bio.gri")
#Env = stack("~/bio.gri")
#Env = getData("worldclim",var="bio",res=2.5)
#Env
#ext = extent(c(-100, - 60, 25, 55))
ext = extent(c(-120, -80, 15, 40))
Env = crop(Env, ext)

#max = maxent(Env, loc) ## dismo maxent function with defaults
## running the model
#?maxent()
#max
#pred = predict(Env, max)
## casting to geo space

## make sure to get rid of NA values
extr = extract(Env, loc) ## gets values from Env that are at loc 
#head(extr)
loc = loc[!is.na(extr[,1]),] ## removes any points in loc with NA values in extr
#nrow(loc)

res = ENMevaluate(occ=loc, env = Env, method='block', 
                  parallel=T, numCores=4, fc=c("L", "LQ", "H"), 
                  RMvalues=seq(0.5,4,0.5), rasterPreds=F)
## this uses your occurences and trimmed env data. it will
## separate the data into 4 blocks, split sample into 4 equal quads
## running in parallel, with 4 cores. 
## the feature classes used are Linear, Linear+Quadratic, and Hunge.
## the regularization values will range from 1/2 to 4 by 1/2s.
## and rasterPreds=F means not running the predict function we ran before
## and thus can't have any AICc

#print(res@results)
#table = res@results
## Mean.AUC is the one to use, 
## full.aUC it's the AUC of the model used for AICc comparisons
## using all the points and no withheld data 
## for omission rate, if small sample size and high confidence
## do mean.ORmin. if don't know, do mean.OR10

## now find the best models 
## minimize omission rate and optimize AUC values
## remever to change ORmin to OR10 depending on sample size 
setsort = res@results[order(res@results[,'Mean.ORmin']),]
setsort2 = setsort[order(setsort[,'Mean.AUC'], decreasing=TRUE),]
top = setsort2[1,]
#top

## CARDINALIS
"  settings features  rm full.AUC  Mean.AUC     Var.AUC Mean.AUC.DIFF
2   LQ_0.5       LQ 0.5   0.9148 0.8975014 0.007018451    0.03151321
Var.AUC.DIFF Mean.OR10   Var.OR10 Mean.ORmin Var.ORmin AICc delta.AICc w.AIC
2  0.004904318 0.1645806 0.02463436      0.022  0.001936   NA         NA    NA
nparam
2     26"
## MESQUITE
"24      H_4        H  4   0.9607 0.9362948 0.008247966    0.02946219
Var.AUC.DIFF Mean.OR10   Var.OR10 Mean.ORmin Var.ORmin AICc delta.AICc w.AIC
24  0.006909553 0.1841774 0.08409337      0.004   6.4e-05   NA         NA    NA
nparam
24     67"

best = which(as.character(res@results[,1]) == as.character(setsort2[1,1]))
pred.raw = predict(Env, res@models[[best]])
writeRaster(pred.raw,filename="PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
png("PredRaw_BestModel_glandulosa.png")
plot(pred.raw, col=viridis::viridis(99))
dev.off()

## thin the data
df = data.frame(occdf)
extr = extract(Env, cbind(df$longitude, df$latitude))
df = df[!is.na(extr[,1]),]

 thin<-thin(loc.data = df, 
            lat.col = "latitude", 
            long.col = "longitude",
            spec.col = "name", 
            thin.par = 10, ## km distance that records need to be separated by
            reps = 10, ## number of times to repeat thinning process
            locs.thinned.list.return = T, 
            write.files = T, 
            max.files = 2, 
            out.dir = "ThinlocsMesquite/",
            #out.dir = "Thinlocs/"
            write.log.file = T)



## import your thinned data
newthin = read.csv('ThinlocsMesquite/thinned_data_thin1.csv')
newthin2 = read.csv("Thinlocs/thinned_data_thin1.csv")

extr = extract(Env, newthin[,2:3]) ## gets values from Env that are at loc 
head(extr)
newthin = newthin[!is.na(extr[,1]),] ## removes any points in loc with NA values in extr
nrow(newthin)

## re-run enm eval on the thinned data
thinres = ENMevaluate(occ=newthin[,2:3], env = Env, 
                      method='block', parallel=T, 
                      numCores=4, fc=c("L", "LQ", "H"),
                      RMvalues=seq(1,4,1), rasterPreds=F)
print(thinres@results)
table2 = thinres@results

png("ThinnedPlot_Mean.AUC_glandulosa.png")
eval.plot(thinres@results, "Mean.AUC")
dev.off()

#and for the thinned model...
setsort = thinres@results[order(thinres@results[,'Mean.ORmin']),]
head(setsort)
setsort2 = setsort[order(setsort[,'Mean.AUC'], decreasing=TRUE),]
head(setsort2)
top = setsort2[1,]
top

## CARDINALIS
"  settings features rm full.AUC  Mean.AUC    Var.AUC Mean.AUC.DIFF Var.AUC.DIFF
9      H_3        H  3   0.9014 0.8484116 0.02874788    0.06941165   0.02755969
Mean.OR10  Var.OR10 Mean.ORmin  Var.ORmin AICc delta.AICc w.AIC nparam
9      0.25 0.1036501 0.05681818 0.01291322   NA         NA    NA     62"

## MESQUITE
" settings features rm full.AUC  Mean.AUC    Var.AUC Mean.AUC.DIFF Var.AUC.DIFF
9      H_3        H  3   0.9014 0.8484116 0.02874788    0.06941165   0.02755969
  Mean.OR10  Var.OR10 Mean.ORmin  Var.ORmin AICc delta.AICc w.AIC nparam
9      0.25 0.1036501 0.05681818 0.01291322   NA         NA    NA     62"

best.thin = which(as.character(thinres@results[,1]) == as.character(setsort2[1,1]))
pred.thin = predict(Env, thinres@models[[best.thin]])
writeRaster(pred.thin,filename="Thinned_PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
png("Thinned_PredRaw_BestModel_glandulosa.png")
plot(pred.thin, col=viridis::viridis(99))
points(newthin[,2:3],col="red")
dev.off()

#print(res@results[best,])
#print(thinres@results[best.thin,])

## now do thresholding
## get the enmeval output
ev.set <- evaluate(newthin[,2:3], thinres@bg.pts, thinres@models[[best.thin]], Env)
th1 = threshold(ev.set) ## omission options
p1.nomit = pred.thin >= th1$no_omission ## the highest thresh where no points are omittec
p1.equal = pred.thin >= th1$equal_sens_spec ## equal sensitivity and specificity according to ROC curve

writeRaster(p1.nomit,filename="Thinned_Thresh_NoOmission_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p1.equal,filename="Thinned_Thresh_EqualSensSpec_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
png("Thinned_Thresh_NoOmission_BestModel_glandulosa.png")
plot(p1.nomit, col=viridis::viridis(99))
dev.off()
png("Thinned_Thresh_EqualSensSpec_BestModel_glandulosa.png")
plot(p1.equal, col=viridis::viridis(99))
dev.off()

## do niche overlap, just between raw and thin niches
nicheOverlap(pred.raw, pred.thin, "D", mask=T)
## for cardinalis sinuatus,0.7759149
## for mesquite, 0.5187202 

## use min convex polygon to do background data
## first generate minimum convex polygon
mcp <- function (xy) { 
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}

MCP.thinlocs = mcp(df[,2:3])
png("MinConvexPolygon_glandulosa.png")
plot(Env[[1]], col=viridis::viridis(99))
plot(MCP.thinlocs, add=T)
dev.off()

## get the background points
envPoly <- rasterToPolygons(Env[[1]], fun=NULL, na.rm=T)

#Get bg
bg.area.thinlocs <- gIntersection(envPoly, MCP.thinlocs)

## turn these into bg points
MCP.raster.thinlocs <- rasterize(bg.area.thinlocs, Env[[1]])
bg.points.thinlocs <- randomPoints(mask = MCP.raster.thinlocs, n = 5000)

## clip just to landmasses
png("MinConvexPolygon_points_glandulosa.png")
plot(Env[[1]], col=viridis::viridis(99))
plot(bg.area.thinlocs, add=T)
dev.off()

## use this to give background points to enmeval fxn, see fxn to do
## use bg.points.thinlocs to build a new best model

## most likely need to change back to newthin and then put bgpoints 
bgpointsres = ENMevaluate(occ=newthin[,2:3],
                          bg.coords=bg.points.thinlocs, env = Env, 
                      method='block', parallel=T, 
                      numCores=4, fc=c("L", "LQ", "H"),
                      RMvalues=seq(1,4,1), rasterPreds=F)
print(bgpointsres@results)
table3 = bgpointsres@results

png("BgPointsPlot_Mean.AUC_glandulosa.png")
eval.plot(bgpointsres@results, "Mean.AUC")
dev.off()

setsort = bgpointsres@results[order(bgpointsres@results[,'Mean.ORmin']),]
head(setsort)
setsort2 = setsort[order(setsort[,'Mean.AUC'], decreasing=TRUE),]
head(setsort2)
top = setsort2[1,]

#sinuatus
"   settings features rm full.AUC  Mean.AUC    Var.AUC Mean.AUC.DIFF
10      L_4        L  4   0.7409 0.6854687 0.05871278     0.1104313
Var.AUC.DIFF Mean.OR10   Var.OR10 Mean.ORmin    Var.ORmin AICc delta.AICc
10   0.02912993 0.2159091 0.04838154 0.01136364 0.0001721763   NA         NA
w.AIC nparam
10    NA      8"
## mesquite
"   settings features rm full.AUC  Mean.AUC    Var.AUC Mean.AUC.DIFF
10      L_4        L  4   0.7409 0.6854687 0.05871278     0.1104313
Var.AUC.DIFF Mean.OR10   Var.OR10 Mean.ORmin    Var.ORmin AICc delta.AICc
10   0.02912993 0.2159091 0.04838154 0.01136364 0.0001721763   NA         NA
w.AIC nparam
10    NA      8"

best.bg = which(as.character(bgpointsres@results[,1]) == as.character(setsort2[1,1]))
pred.bg = predict(Env, bgpointsres@models[[best.bg]])
writeRaster(pred.bg,filename="BgPointsPlot_PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
png("BgPointsPlot_PredRaw_BestModel_glandulosa.png")
plot(pred.bg, col=viridis::viridis(99))
#points(bg.points.thinlocs,col="red")
dev.off()

print(bgpointsres@results[best.bg,])

ev.set <- evaluate(newthin[,2:3], bgpointsres@bg.pts, 
                   bgpointsres@models[[best.bg]], Env)
th2 = threshold(ev.set) ## omission options
p2.nomit = pred.bg >= th2$no_omission ## the highest thresh where no points are omittec
p2.equal = pred.bg >= th2$equal_sens_spec ## equal sensitivity and specificity according to ROC curve

writeRaster(p2.nomit,filename="BgPoints_Thresh_NoOmission_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p2.equal,filename="BgPoints_Thresh_EqualSensSpec_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
png("BgPoints_Thresh_NoOmission_BestModel_glandulosa.png")
plot(p2.nomit, col=viridis::viridis(99))
dev.off()
png("BgPoints_Thresh_EqualSensSpec_BestModel_glandulosa.png")
plot(p2.equal, col=viridis::viridis(99))
dev.off()

## weird results
## TO DO ##
## if you add polygon is it doing better in region it occurs or not
## try buffered sampling for this spp, 500 km around each point? 

## try using a differnet dataset instead of worldclim
## change Env to envirem?
enirem = stack('/data/spbio/climgrids/envirem.gri')
enirem
#ext = extent(c(-100, - 60, 25, 55))
ext = extent(c(-120, - 80, 15, 40))
enirem = crop(enirem, ext)

extr = extract(Env, newthin[,2:3])

enviremres = ENMevaluate(occ=newthin[,2:3], env = enirem, 
                          method='block', parallel=T, 
                          numCores=4, fc=c("L", "LQ", "H"),
                          RMvalues=seq(1,4,1), rasterPreds=F)
print(enviremres@results)
table4 = enviremres@results

png("EnviremPlot_Mean.AUC_glandulosa.png")
eval.plot(enviremres@results, "Mean.AUC")
dev.off()

setsort = enviremres@results[order(enviremres@results[,'Mean.ORmin']),]
head(setsort)
setsort2 = setsort[order(setsort[,'Mean.AUC'], decreasing=TRUE),]
head(setsort2)
top = setsort2[1,]
##settings features rm full.AUC  Mean.AUC    Var.AUC Mean.AUC.DIFF
##11     LQ_4       LQ  4   0.8929 0.8536515 0.02662025    0.06999511
##Var.AUC.DIFF Mean.OR10   Var.OR10 Mean.ORmin   Var.ORmin AICc delta.AICc
##11   0.02002682 0.2580603 0.06241952 0.04637949 0.006137704   NA         NA
##w.AIC nparam
##11    NA     12
best.env = which(as.character(enviremres@results[,1]) == as.character(setsort2[1,1]))
pred.env = predict(enirem, enviremres@models[[best.env]])
writeRaster(pred.env,filename="EnviremPlot_PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
png("EnviremPlot_PredRaw_BestModel_glandulosa.png")
plot(pred.env, col=viridis::viridis(99))
#points(bg.points.thinlocs,col="red")
dev.off()

print(enviremres@results[best.env,])

ev.set <- evaluate(newthin[,2:3], enviremres@bg.pts, 
                   enviremres@models[[best.env]], enirem)
th3 = threshold(ev.set) ## omission options
p3.nomit = pred.env >= th3$no_omission ## the highest thresh where no points are omittec
p3.equal = pred.env >= th3$equal_sens_spec ## equal sensitivity and specificity according to ROC curve

writeRaster(p3.nomit,filename="Envirem_Thresh_NoOmission_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p3.equal,filename="Envirem_Thresh_EqualSensSpec_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)

png("Envirem_Thresh_NoOmission_BestModel_glandulosa.png")
plot(p3.nomit, col=viridis::viridis(99))
dev.off()
png("Envirem_Thresh_EqualSensSpec_BestModel_glandulosa.png")
plot(p3.equal, col=viridis::viridis(99))
dev.off()

## and lastly do paleoclimate models
#lgm = stack('/data/spbio/climgrids/miroc/lgm/bio.gri')
#lgm.pred.thin = predict(lgm,  thinres@models[[best.thin]]) 

## first do microc
lgm = stack('/data/spbio/climgrids/miroc/lgm_wc.gri')
## need to change miroc to mpi and ccsm4 as needed
names(lgm) = names(Env) #Make names match 
lgm = crop(lgm, extent(Env))
lgm.pred.thin = predict(lgm,  thinres@models[[best.thin]]) 
lgm.pred.thin


png("LGM_Microc_PredRaw_BestModel_glandulosa.png")
plot(lgm.pred.thin, col=viridis::viridis(99))
#points(newthin[,2:3],col="red")
dev.off()

ev.set <- evaluate(newthin[,2:3], thinres@bg.pts, 
                   thinres@models[[best.env]], lgm)
th4 = threshold(ev.set) ## omission options
p4.nomit = pred.env >= th4$no_omission ## the highest thresh where no points are omittec
p4.equal = pred.env >= th4$equal_sens_spec ## equal sensitivity and specificity according to ROC curve

writeRaster(lgm.pred.thin,filename="LGM_Microc_PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p4.nomit,filename="LGM_Miroc_Thresh_NoOmission_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p4.equal,filename="LGM_Miroc_Thresh_EqualSensSpec_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)

png("LGM_Miroc_Thresh_NoOmission_BestModel_glandulosa.png")
plot(p4.nomit, col=viridis::viridis(99))
dev.off()
png("LGM_Miroc_Thresh_EqualSensSpec_BestModel_glandulosa.png")
plot(p4.equal, col=viridis::viridis(99))
dev.off()

## now mpi
lgm = stack('/data/spbio/climgrids/mpi/lgm_wc.gri')
names(lgm) = names(Env) #Make names match 
lgm = crop(lgm, extent(Env))
lgm.pred.thin = predict(lgm,  thinres@models[[best.thin]]) 
lgm.pred.thin

png("LGM_Mpi_PredRaw_BestModel_glandulosa.png")
plot(lgm.pred.thin, col=viridis::viridis(99))
#points(newthin[,2:3],col="red")
dev.off()

ev.set <- evaluate(newthin[,2:3], thinres@bg.pts, 
                   thinres@models[[best.env]], lgm)
th5 = threshold(ev.set) ## omission options
p5.nomit = pred.env >= th5$no_omission ## the highest thresh where no points are omittec
p5.equal = pred.env >= th5$equal_sens_spec ## equal sensitivity and specificity according to ROC curve

writeRaster(lgm.pred.thin,filename="LGM_Mpi_PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p5.nomit,filename="LGM_Mpi_Thresh_NoOmission_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p5.equal,filename="LGM_Mpi_Thresh_EqualSensSpec_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)


png("LGM_Mpi_Thresh_NoOmission_BestModel_glandulosa.png")
plot(p5.nomit, col=viridis::viridis(99))
dev.off()
png("LGM_Mpi_Thresh_EqualSensSpec_BestModel_glandulosa.png")
plot(p5.equal, col=viridis::viridis(99))
dev.off()

## now ccsm4
lgm = stack('/data/spbio/climgrids/ccsm4/lgm_wc.gri')
names(lgm) = names(Env) #Make names match 
lgm = crop(lgm, extent(Env))
lgm.pred.thin = predict(lgm,  thinres@models[[best.thin]]) 
lgm.pred.thin

png("LGM_Ccsm4_PredRaw_BestModel_glandulosa.png")
plot(lgm.pred.thin, col=viridis::viridis(99))
#points(newthin[,2:3],col="red")
dev.off()

ev.set <- evaluate(newthin[,2:3], thinres@bg.pts, 
                   thinres@models[[best.env]], lgm)
th6 = threshold(ev.set) ## omission options
p6.nomit = pred.env >= th6$no_omission ## the highest thresh where no points are omittec
p6.equal = pred.env >= th6$equal_sens_spec ## equal sensitivity and specificity according to ROC curve

writeRaster(lgm.pred.thin,filename="LGM_Ccsm4_PredRaw_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p6.nomit,filename="LGM_Ccsm4_Thresh_NoOmission_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)
writeRaster(p6.equal,filename="LGM_Ccsm4_Thresh_EqualSensSpec_BestModel_glandulosa.asc",
            format="ascii",overwrite=T)

png("LGM_Ccsm4_Thresh_NoOmission_BestModel_glandulosa.png")
plot(p6.nomit, col=viridis::viridis(99))
dev.off()
png("LGM_Ccsm4_Thresh_EqualSensSpec_BestModel_glandulosa.png")
plot(p6.equal, col=viridis::viridis(99))
dev.off()

##Homework 1
## For tomorrow:
##Pick two (or more) new species from a single genus (or co-occurring/interacting species) 
##and use ENMeval to find the best model for each.
##Make maps of best model output for each as continuous and binary ranges.
##Record your "best" model parameters and be able to discuss all of your 
##modeling choices (predictor variables, threshold, thinning, model evaluation, etc.)
##Put this in a powerpoint slide and email to Pete by 9:30 AM tomorrow morning




## create a dist model for a spp of interest to you using enmeval, make a few slides
##Pick two (or more) new species from a single genus and use ENMeval to find the best model for each.
##Make maps of best model output for each as continuous and binary ranges.
## do model parameters, why you picked stuff, etc. make sure thinned. 
## rerun all this code once done to make sure going well
##Project into LGM for all 3 GCMs. (2 sp. x 3 GCMs)
## Calculate range overlap between each LGM projection and the modern as a proxy for geographic range shift.
##Calculate range overlap between each LGM projection and the modern as a proxy for geographic range shift.





## save
save.image("out.RData")








