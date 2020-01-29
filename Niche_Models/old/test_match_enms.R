layer1 = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/USA_ONLY/Toxostoma/crissale/AddedPredThin_BestModel_Toxostoma crissale_addedLayers_USA.asc_-4.63216781616211_-17.7025375366211.rescaled.asc"
layer2 = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/WORLDCLIM/Toxostoma/crissale/AddedPredThin_BestModel_Toxostoma crissale_addedLayers_worldclim.asc_-4.94697666168213_-46.4806938171387.rescaled.asc"

l1 = raster::raster(layer1)
l2 = raster::raster(layer2)

l3=crop(l1,extent(l2))
l4 = resample(l3,l2)

v3 = scales::rescale(values(l2))
v4 = scales::rescale(values(l4))

v5 = cbind(v3,v4)
v6 = v5[complete.cases(v5),]

v6[,1] = scales::rescale(v6[,1])
v6[,2] = scales::rescale(v6[,2])
colnames(v6) = c("Deserts","Continent")
v6 = as.data.frame(v6)

cor(v6,use="pairwise.complete.obs")^2
png("CRISSALE-MATCH-ENMS.png")
plot(v6$Deserts,v6$Continent)
abline(lm(v6$Continent~v6$Deserts),col="red")
dev.off()
