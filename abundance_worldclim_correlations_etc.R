worldclim_enm = "/Users/kprovost/Downloads/worldclimasci/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_worldclim.asc_-6.46640157699585_-69.5545196533203.rescaled.asc"
continent_enm = "/Users/kprovost/Downloads/continent asci/AddedPredThin_BestModel_Amphispiza bilineata_addedLAyers.asc_-3.78692746162415_-58.8649635314941.rescaled.asc"
abundance_asc = "/Users/kprovost/Downloads/clipped abun/BILINEATA_Idw_abundance_clipped.asc"
raw_env_tif = "/Users/kprovost/Downloads/ENMS_multilayer4.tif"

worldclim = raster(worldclim_enm)
continent = raster(continent_enm)
abundance = raster(abundance_asc)
allrawenv = stack(raw_env_tif)

res(worldclim)
res(continent)
res(abundance)
res(allrawenv)

wc_pt = rasterToPoints(worldclim, spatial = F)
ct_pt = rasterToPoints(continent, spatial = F)
ab_pt = rasterToPoints(abundance, spatial = F)

colnames(wc_pt) = c("x","y","wc")
colnames(ct_pt) = c("x","y","ct")
colnames(ab_pt) = c("x","y","ab")

merged = merge(wc_pt,ct_pt,by=c("x","y"))
merged = merge(merged,ab_pt,by=c("x","y"))
dim(merged)
en_pt = rasterToPoints(allrawenv, spatial = F)

#colnames(en_pt) = c("x","y","en")
cor(merged[,c("wc","ab")],use="pairwise.complete.obs")
png("wc_ab_bil.png")
plot(merged$wc,merged$ab)
dev.off()

