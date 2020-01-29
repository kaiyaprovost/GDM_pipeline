library(raster)

abd_path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/abundance_asciis/clipped/"
enm_path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/3.enms/USA_ONLY/"

abd_list = c("BELLII_Idw_abundance_clipped.asc",
             "CURVIROSTRE_Idw_abundance_clipped.asc",
             "CRISSALE_Idw_abundance_clipped.asc",
             "MELANURA_Idw_abundance_clipped.asc",
             "NITENS_Idw_abundance_clipped.asc",
             "FUSCA_Idw_abundance_clipped.asc",
             "SINUATUS_Idw_abundance_clipped.asc",
             "BRUNNEICAPILLUS_Idw_abundance_clipped.asc",
             "FLAVICEPS_Idw_abundance_clipped.asc",
             "BILINEATA_Idw_abundance_clipped.asc")

enm_list = c("Vireo/AddedPredThin_BestModel_Vireo bellii_addedLayers_USA.asc_-4.56468915939331_-18.812442779541.rescaled.asc",
             "Toxostoma/cur_full/AddedPredThin_BestModel_Toxostoma curvirostre_addedLayers_USA.asc_-5.62186193466187_-12.8309497833252.rescaled.asc",
             "Toxostoma/crissale/AddedPredThin_BestModel_Toxostoma crissale_addedLayers_USA.asc_-4.63216781616211_-17.7025375366211.rescaled.asc",
             "Polioptila/AddedPredThin_BestModel_Polioptila melanura_addedLayers_USA.asc_-3.93294405937195_-21.4272289276123.rescaled.asc",
             "Phainopepla/AddedPredThin_BestModel_Phainopepla nitens_addedLayers_USA.asc_-3.93656158447266_-22.6519756317139.rescaled.asc",
             "Melozone/AddedPredThin_BestModel_Melozone fusca_addedLayers_USA.asc_-5.73480558395386_-13.4299039840698.rescaled.asc",
             "Cardinalis/AddedPredThin_BestModel_Cardinalis sinuatus_addedLayers_USA.asc_-5.24065685272217_-19.4379730224609.rescaled.asc",
             "Campylorhynchus/AddedPredThin_BestModel_Campylorhynchus brunneicapillus_addedLayers_USA.asc_-5.1320424079895_-18.1571197509766.rescaled.asc",
             "Auriparus/AddedPredThin_BestModel_Auriparus flaviceps_addedLayers_USA.asc_-4.34592342376709_-18.2621974945068.rescaled.asc",
             "Amphispiza/AddedPredThin_BestModel_Amphispiza bilineata_addedLayers_USA.asc_-5.90643692016602_-13.0102529525757.rescaled.asc")

for (i in 1:10) {
print(i)
abd = raster(paste(abd_path,abd_list[i],sep=""))
enm = raster(paste(enm_path,enm_list[i],sep=""))

print(abd)
print(enm)

ext_a = extent(abd)
ext_e = extent(enm)

abd_c = crop(abd,ext_e)

enm_c = aggregate(enm,20)
#abd_cr = disaggregate(abd_c,20)

abd_val = values(abd_c)
enm_val = values(enm_c)

val = data.frame(abd_val,enm_val)
val = unique(val)
val = val[complete.cases(val),]
names(val) = c("abd","enm")

nrow(val)
cor(val)

subone = val[val$abd < 1,]
moreone = val[val$abd >= 1,]

#subonequan = quantile(subone$enm,probs=c(0.25,0.50,0.75))
moreonequan = quantile(moreone$enm,probs=c(0.25,0.50,0.75))

always = max(subone$enm)

never = min(moreone$enm)

png(paste(abd_path,abd_list[i],".HIST.png",sep=""))
hist(val$enm,main=abd_list[i],xlab="Environmental Value")
abline(v=always,col="red")
abline(v=never,col="blue")
dev.off()

png(paste(abd_path,abd_list[i],".RASTER.png",sep=""))
breakpoints <- c(0,never,moreonequan,always,1)
colors <- c("black","red","pink","yellow","lightblue","blue")
plot(enm,breaks=breakpoints,col=colors,main=abd_list[i])
legend("bottomleft",legend=c("Never","<25%","25-50%","50-75%","75+%","Always"),
       fill=colors)
dev.off()

png(paste(abd_path,abd_list[i],".png",sep=""))
plot(val$abd,val$enm,col=rgb(0,0,0,0.1),
     main=abd_list[i],xlab="Abundance",
     ylab="Habitat Suitability")
mod = lm(val$enm~sqrt(sqrt(val$abd)))
summary(mod)
rsq = summary(mod)$adj.r.squared
newx = (val$abd)
newy = predict(mod,data.frame(abd=newx))
points(newx,newy, col=rgb(1,0,0,0.1))
legend("bottomright",legend=c(paste("Rsq: ",round(rsq,2),sep=""),
                              paste("Always: ",round(always,2),sep=""),
                              paste("Never: ",round(never,2),sep="")),
       bty="n")
abline(h=always,col="blue",lty=2,lwd=2)
abline(h=never,col="blue",lty=2,lwd=2)

dev.off()

}

