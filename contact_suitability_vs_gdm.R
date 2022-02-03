species=c("BELLII","BILINEATA","BRUNNEICAPILLUS","CRISSALE","CURVIROSTRE","FLAVICEPS","FUSCA","MELANURA","NITENS","SINUATUS")
allmean=c(0.65,0.54,0.65,0.52,0.57,0.57,0.55,0.58,0.48,0.71)
allsd=c(0.15,0.12,0.14,0.17,0.14,0.15,0.13,0.15,0.18,0.15)
IBA=c(0.06,0.10,0.58,0.25,0.23,0.24,0.14,0.12,0.06,0.04)
IBB=c(0.13,0.16,0.06,0.23,0.38,0.07,0.34,0.20,0.04,0.56)
IBD=c(0.06,0.16,0.16,0.17,0.13,0.26,0.12,0.16,0.04,0.18)
IBE=c(0.45,0.51,0.10,0.19,0.21,0.31,0.28,0.49,0.27,0.18)
IBH=c(0.11,0.02,0.08,0.17,0.00,0.11,0.10,0.04,0.49,0.04)
MIXED=c(0.19,0.06,0.02,0.00,0.06,0.00,0.02,0.00,0.10,0.00)

means=c(0.59,0.60,0.59,0.58,0.54,0.58)
sds=c(0.15,0.14,0.15,0.15,0.16,0.15)
models=c("IBA","IBB","IBD","IBE","IBH","MIXED")
par(mfrow=c(1,1))
plot(means,ylim=c(0,1),type="b")
points(means-sds,col="red",type="b")
points(means+sds,col="red",type="b")

df = read.table("/Users/kprovost/Documents/GitHub/bioacoustics/tweetynet/EV/test_contact_suitability_vs_gdms.txt",
                sep="\t",header=T)

pdf("~/Dropbox (AMNH)/suitability_mean_sd_contact_zone.pdf")

par(mfrow=c(2,3))
boxplot(df$contact_suitability_mean~df$MOSTA.1MODEL1)
mod=aov(df$contact_suitability_mean~df$MOSTA.1MODEL1)
summary(mod)
TukeyHSD(mod)
mod=aov(df$contact_suitability_mean[df$SPECIES!="NITENS"]~df$MOSTA.1MODEL1[df$SPECIES!="NITENS"])
summary(mod)
TukeyHSD(mod)
boxplot(df$contact_suitability_mean~df$MOSTA.1MODEL2)
boxplot(df$contact_suitability_mean[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
boxplot(df$contact_suitability_sd~df$MOSTA.1MODEL1)
mod=aov(df$contact_suitability_sd~df$MOSTA.1MODEL1)
summary(mod)
TukeyHSD(mod)
mod=aov(df$contact_suitability_sd[df$SPECIES!="NITENS"]~df$MOSTA.1MODEL1[df$SPECIES!="NITENS"])
summary(mod)
TukeyHSD(mod)
boxplot(df$contact_suitability_sd~df$MOSTA.1MODEL2)
boxplot(df$contact_suitability_sd[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])

par(mfrow=c(1,1))
plot(1:10,1:10,col=1:10,pch=1:10,cex=3)
legend("topleft",legend=species,col=1:10,pch=1:10,cex=2)

par(mfrow=c(2,1))
plot(IBA,allmean,col=1:10,pch=1:10,main="IBA")
mod=lm(allmean~IBA)
summary(mod)
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBA,allsd,col=1:10,pch=1:10)
mod2=lm(allsd~IBA)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(IBB,allmean,col=1:10,pch=1:10,main="IBB")
mod=lm(allmean~IBB)
summary(mod)
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBB,allsd,col=1:10,pch=1:10)
mod2=lm(allsd~IBB)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(IBB[1:9],allmean[1:9],col=1:9,pch=1:9,main="IBB-nosin")
mod=lm(allmean[1:9]~IBB[1:9])
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBB[1:9],allsd[1:9],col=1:9,pch=1:9)
mod2=lm(allsd[1:9]~IBB[1:9])
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(IBD,allmean,col=1:10,pch=1:10,main="IBD")
mod=lm(allmean~IBD)
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBD,allsd,col=1:10,pch=1:10)
mod2=lm(allsd~IBD)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(IBE,allmean,col=1:10,pch=1:10,main="IBE")
mod=lm(allmean~IBE)
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBE,allsd,col=1:10,pch=1:10)
mod2=lm(allsd~IBE)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(IBH,allmean,col=1:10,pch=1:10,main="IBH")
mod=lm(allmean~IBH)
summary(mod)
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBH,allsd,col=1:10,pch=1:10)
mod2=lm(allsd~IBH)
summary(mod2)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(IBH[c(1:8,10)],allmean[c(1:8,10)],col=c(1:8,10),pch=c(1:8,10),main="IBH-nonit")
mod=lm(allmean[c(1:8,10)]~IBH[c(1:8,10)])
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(IBH[c(1:8,10)],allsd[c(1:8,10)],col=c(1:8,10),pch=c(1:8,10))
mod2=lm(allsd[c(1:8,10)]~IBH[c(1:8,10)])
summary(mod2)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

par(mfrow=c(2,1))
plot(MIXED,allmean,col=1:10,pch=1:10,main="MIXED")
mod=lm(allmean~MIXED)
abline(mod,col="black",lty=3)
mtext(round(summary(mod)$adj.r.squared,2))
plot(MIXED,allsd,col=1:10,pch=1:10)
mod2=lm(allsd~MIXED)
abline(mod2,col="red",lty=3)
mtext(round(summary(mod2)$adj.r.squared,2))

dev.off()


corr = read.table("/Users/kprovost/Dropbox (AMNH)/correlation_Between_reconmbination_gdmscale.txt",
                  sep="\t",header=T)
corr=corr[,c("DATASET","BELLII","FLAVICEPS","NITENS","MELANURA","BRUNNEICAPILLUS","CRISSALE","CURVIROSTRE","BILINEATA","FUSCA","SINUATUS")]
calccor=cor(corr[,2:11],use="pairwise.complete.obs")
corrplot::corrplot(calccor,method="color",
                   col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(200))




rec = data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_fix7.temp",
                        data.table = F)
rec = rec[,c("species","chrpos","weighted_recomb")]
rec=unique(rec)

bel = rec[rec$species=="bel",c("chrpos","weighted_recomb")]; colnames(bel) = c("chrpos","bel")
bil = rec[rec$species=="bil",c("chrpos","weighted_recomb")]; colnames(bil) = c("chrpos","bil")
bru = rec[rec$species=="bru",c("chrpos","weighted_recomb")]; colnames(bru) = c("chrpos","bru")
cri = rec[rec$species=="cri",c("chrpos","weighted_recomb")]; colnames(cri) = c("chrpos","cri")
cur = rec[rec$species=="cur",c("chrpos","weighted_recomb")]; colnames(cur) = c("chrpos","cur")
fla = rec[rec$species=="fla",c("chrpos","weighted_recomb")]; colnames(fla) = c("chrpos","fla")
fus = rec[rec$species=="fus",c("chrpos","weighted_recomb")]; colnames(fus) = c("chrpos","fus")
mel = rec[rec$species=="mel",c("chrpos","weighted_recomb")]; colnames(mel) = c("chrpos","mel")
nit = rec[rec$species=="nit",c("chrpos","weighted_recomb")]; colnames(nit) = c("chrpos","nit")
sin = rec[rec$species=="sin",c("chrpos","weighted_recomb")]; colnames(sin) = c("chrpos","sin")

d1 = merge(bel,bil,all=T)
d2 = merge(bru,cri,all=T)
d3 = merge(cur,fla,all=T)
d4 = merge(fus,mel,all=T)
d5 = merge(nit,sin,all=T)

d1 = merge(d1,d2,all=T)
d2 = merge(d3,d4,all=T)
d3 = d5

d1 = merge(d1,d2,all=T)
d2 = d5

d1 = merge(d1,d2,all=T)
mergecor = cor(d1[,2:11],use="pairwise.complete.obs")

corrplot::corrplot(mergecor)



