setwd("~")

df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate_gdm_results_useold.csv",
                sep=",",
                header=T,fill=T,
                stringsAsFactors = F,
                skip = 0)

png("summary_stats_by_species_gdm.png",height=4.5,width=6,units = "in",
    res=300)
par(mfrow=c(2,2),
    cex.axis=0.75,cex.lab=0.75,
    mar=c(0.5,4,0.5,0))
cols=rainbow(10)
plotrix::gap.boxplot(df$MEAN_FST~substr(df$SPECIES,1,3),
                     gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                     axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",
                     xaxt="n",col=cols)
title(ylab="Mean Fst")
axis(2,labels=c(seq(0,0.2,0.05),0.57),
     at=c(seq(0,0.2,0.05),0.235),tick=T)
boxplot(df$MEAN_DXY~substr(df$SPECIES,1,3),ylab="Mean Dxy",
        xaxt="n",col=cols)
boxplot(df$PERCENT_MISSING~substr(df$SPECIES,1,3),
        ylab="Percent Missing",xaxt="n",col=cols)
boxplot((df$MEAN_RECOMB/(1e-9))~substr(df$SPECIES,1,3),
        ylab="Mean Recombination Rate (x 1^-9)",
        xaxt="n",col=cols)
dev.off()



numonly = df[,c("MEAN_FST","MEAN_DXY","MEAN_RECOMB","PERCENT_MISSING")]
numonly_nofstout = numonly[numonly$MEAN_FST<=0.5,]
png("numonly_nofstout.png")
plot(numonly_nofstout)
dev.off()

cor(numonly,use="pairwise.complete.obs")
cor(numonly_nofstout,use="pairwise.complete.obs")

png("numonly_plus_numonly_nofstout_corrplots.png")
par(mfrow=c(1,2))
corrplot::corrplot(cor(numonly,use="pairwise.complete.obs"),
                   method = "number",diag=T,type="upper")
corrplot::corrplot(cor(numonly_nofstout,use="pairwise.complete.obs"),
                   method="number",diag=T,type="upper")
dev.off()

png("numonly_correlations.png")
par(mfrow=c(2,3),mar=c(4,4,0,0))
plot(numonly$PERCENT_MISSING,numonly$MEAN_FST)
abline(lm(numonly$MEAN_FST~numonly$PERCENT_MISSING),col="red")
summary(lm(numonly$MEAN_FST~numonly$PERCENT_MISSING))$adj.r.squared

plot(numonly$PERCENT_MISSING,numonly$MEAN_DXY)
abline(lm(numonly$MEAN_DXY~numonly$PERCENT_MISSING),col="red")
summary(lm(numonly$MEAN_DXY~numonly$PERCENT_MISSING))$adj.r.squared

plot(numonly$PERCENT_MISSING,numonly$MEAN_RECOMB)
abline(lm(numonly$MEAN_RECOMB~numonly$PERCENT_MISSING),col="red")
summary(lm(numonly$MEAN_RECOMB~numonly$PERCENT_MISSING))$adj.r.squared

plot(numonly$MEAN_DXY,numonly$MEAN_FST)
abline(lm(numonly$MEAN_FST~numonly$MEAN_DXY),col="red")
summary(lm(numonly$MEAN_FST~numonly$MEAN_DXY))$adj.r.squared

plot(numonly$MEAN_RECOMB,numonly$MEAN_FST)
abline(lm(numonly$MEAN_FST~numonly$MEAN_RECOMB),col="red")
summary(lm(numonly$MEAN_FST~numonly$MEAN_RECOMB))$adj.r.squared

plot(numonly$MEAN_RECOMB,numonly$MEAN_DXY)
abline(lm(numonly$MEAN_DXY~numonly$MEAN_RECOMB),col="red")
summary(lm(numonly$MEAN_DXY~numonly$MEAN_RECOMB))$adj.r.squared
dev.off()


png("numonly_nofstout_correlations.png")
par(mfrow=c(2,3),mar=c(4,4,0,0))
plot(numonly_nofstout$PERCENT_MISSING,numonly_nofstout$MEAN_FST)
abline(lm(numonly_nofstout$MEAN_FST~numonly_nofstout$PERCENT_MISSING),col="red")
summary(lm(numonly_nofstout$MEAN_FST~numonly_nofstout$PERCENT_MISSING))$adj.r.squared

plot(numonly_nofstout$PERCENT_MISSING,numonly_nofstout$MEAN_DXY)
abline(lm(numonly_nofstout$MEAN_DXY~numonly_nofstout$PERCENT_MISSING),col="red")
summary(lm(numonly_nofstout$MEAN_DXY~numonly_nofstout$PERCENT_MISSING))$adj.r.squared

plot(numonly_nofstout$PERCENT_MISSING,numonly_nofstout$MEAN_RECOMB)
abline(lm(numonly_nofstout$MEAN_RECOMB~numonly_nofstout$PERCENT_MISSING),col="red")
summary(lm(numonly_nofstout$MEAN_RECOMB~numonly_nofstout$PERCENT_MISSING))$adj.r.squared

plot(numonly_nofstout$MEAN_DXY,numonly_nofstout$MEAN_FST)
abline(lm(numonly_nofstout$MEAN_FST~numonly_nofstout$MEAN_DXY),col="red")
summary(lm(numonly_nofstout$MEAN_FST~numonly_nofstout$MEAN_DXY))$adj.r.squared

plot(numonly_nofstout$MEAN_RECOMB,numonly_nofstout$MEAN_FST)
abline(lm(numonly_nofstout$MEAN_FST~numonly_nofstout$MEAN_RECOMB),col="red")
summary(lm(numonly_nofstout$MEAN_FST~numonly_nofstout$MEAN_RECOMB))$adj.r.squared

plot(numonly_nofstout$MEAN_RECOMB,numonly_nofstout$MEAN_DXY)
abline(lm(numonly_nofstout$MEAN_DXY~numonly_nofstout$MEAN_RECOMB),col="red")
summary(lm(numonly_nofstout$MEAN_DXY~numonly_nofstout$MEAN_RECOMB))$adj.r.squared
dev.off()




png("nmissing_correlations.png")
par(mfrow=c(3,1))
plot(df$MEAN_FST,df$PERCENT_MISSING)
plot(df$MEAN_DXY,df$PERCENT_MISSING)
plot(df$MEAN_RECOMB,df$PERCENT_MISSING)
dev.off()

cor(df$MEAN_FST,df$PERCENT_MISSING,use="pairwise.complete.obs") # -0.337
cor(df$MEAN_DXY,df$PERCENT_MISSING,use="pairwise.complete.obs") # 0.267
cor(df$MEAN_RECOMB,df$PERCENT_MISSING,use="pairwise.complete.obs") # 0.08

library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
png("boxplots.png")
par(mar=c(4,4,0.1,0.1),
    mfrow=c(1,3))
boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""]/(1e-9)~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)

boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL1!=""]/(1e-9)~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""]/(1e-9)~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
dev.off()

model=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
summary(model) 
## sig uni, mixed different
## not significant biv
## almost sig multim mixed-inh almost sig
TukeyHSD(model)

agg=aggregate(df$PERCENT_MISSING~df$SPECIES,FUN=function(x){sd(x,na.rm=T)})


png("four_panel_figure_chapter_3_bivariate.png",height=4,width=6,units = "in",
    res=300)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
## four panels to see if same
par(mfrow=c(2,2),cex.axis=1)
par(mar=c(0.1,4,0.3,0.1))
#boxplot(df$MEAN_FST[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
#        col=cols,xaxt="n",
#        ylab="Mean Fst",
#        xlab="Best Model",las=2)
plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
                     gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                     axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                     xaxt="n",col=cols)
#legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
#       col=cols,fill=cols)
title(ylab="Mean Fst")
axis(2,labels=c(seq(0,0.2,0.05),0.57),
     at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
boxplot(df$MEAN_DXY[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Mean Dxy",
        xlab="Best Model",las=1)
boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""]*100~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Percent Missing",
        xlab="Best Model",las=1)
boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""]/(1e-9)~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Mean Recombination Rate",
        xlab="Best Model",las=2)
dev.off()
## sig uni, mixed different
## not significant biv
## almost sig multim mixed-inh almost sig


png("four_panel_figure_chapter_3_univariate.png",height=4,width=6,units = "in",
    res=300)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,2,3,7,4)],"grey")
## four panels to see if same
par(mfrow=c(2,2))
par(mar=c(0.1,4,0.3,0.1))
#boxplot(df$MEAN_FST[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
#        col=cols,xaxt="n",
#        ylab="Mean Fst",
#        xlab="Best Model",las=2)
plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
                     gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                     axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                     xaxt="n",col=cols)
#legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
#       col=cols,fill=cols)
title(ylab="Mean Fst")
axis(2,labels=c(seq(0,0.2,0.05),0.57),
     at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
boxplot(df$MEAN_DXY[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,xaxt="n",
        ylab="Mean Dxy",
        xlab="Best Model",las=2)
boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL1!=""]*100~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,xaxt="n",
        ylab="Percent Missing",
        xlab="Best Model",las=2)
boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL1!=""]/(1e-9)~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,xaxt="n",
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
dev.off()




png("four_panel_figure_chapter_3_multivariate.png",height=4,width=6,units = "in",
    res=300)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,7,4)],"grey")
## four panels to see if same
par(mfrow=c(2,2))
par(mar=c(0.1,4,0.3,0.1))
#boxplot(df$MEAN_FST[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
#        col=cols,xaxt="n",
#        ylab="Mean Fst",
#        xlab="Best Model",las=2)
plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
                     gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                     axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                     xaxt="n",col=cols)
#legend("topright",legend=c("IBA","IBE","IBH","MIX"),
#       col=cols,fill=cols)
title(ylab="Mean Fst")
axis(2,labels=c(seq(0,0.2,0.05),0.57),
     at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
boxplot(df$MEAN_DXY[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,xaxt="n",
        ylab="Mean Dxy",
        xlab="Best Model",las=2)
boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""]*100~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,xaxt="n",
        ylab="Percent Missing",
        xlab="Best Model",las=2)
boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""]/(1e-9)~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,xaxt="n",
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
dev.off()



#df = df[df$SPECIES!="NITENS",]

## aov model -- BIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])

model1s=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""] + df$SPECIES[df$MOSTA.1MODEL2!=""])
model2s=aov(df$MEAN_FST[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""]+ df$SPECIES[df$MOSTA.1MODEL2!=""])
model3s=aov(df$MEAN_DXY[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""]+ df$SPECIES[df$MOSTA.1MODEL2!=""])
model4s=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""]+ df$SPECIES[df$MOSTA.1MODEL2!=""])
model5s=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""]+ df$SPECIES[df$MOSTA.1MODEL2!=""])


summary(model1) # sig 0.0498
summary(model2) # not
summary(model3) # sig 0.0245
summary(model4) # sig 1.96e-9 -- still sig without nitens
summary(model5) # not sig 0.115

summary(model1s) # sig 0.0498
summary(model2s) # not
summary(model3s) # sig 0.0245
summary(model4s) # sig 1.96e-9 -- still sig without nitens
summary(model5s) # not sig 0.115


TukeyHSD(model1) # nearly ibe-iba 0.0565688
TukeyHSD(model2) # none
TukeyHSD(model3) # sig ibd-iba 0.0340591
TukeyHSD(model4) # ibd-iba 0.0045044 mix-iba 0 mix-ibd 0.0008079 ibh-ibe 0.0042963 mix-ibe 0 mix-ibh 0.0127600
## ibh comparisons not sig without nitens 
TukeyHSD(model5) # ibd-iba 0.036m ibe-ibd 0.066, ibh-ibd = 0.027

TukeyHSD(model1s) # nearly ibe-iba 0.0565688
TukeyHSD(model2s) # none
TukeyHSD(model3s) # sig ibd-iba 0.0340591
TukeyHSD(model4s) # ibd-iba 0.0045044 mix-iba 0 mix-ibd 0.0008079 ibh-ibe 0.0042963 mix-ibe 0 mix-ibh 0.0127600
## ibh comparisons not sig without nitens 
TukeyHSD(model5s) # ibd-iba 0.036m ibe-ibd 0.066, ibh-ibd = 0.027


## aov model -- UNIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])

model1s=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""] + df$SPECIES[df$MOSTA.1MODEL1!=""])
model2s=aov(df$MEAN_FST[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""] + df$SPECIES[df$MOSTA.1MODEL1!=""])
model3s=aov(df$MEAN_DXY[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""] + df$SPECIES[df$MOSTA.1MODEL1!=""])
model4s=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""] + df$SPECIES[df$MOSTA.1MODEL1!=""])
model5s=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""] + df$SPECIES[df$MOSTA.1MODEL1!=""])


summary(model1) # sig 0.0455 -- 30 sep 2021: 0.072
summary(model2) # not 0.846 -- 30 sep 2021: 0.0367
summary(model3) # not 0.527 -- 30 sep 2021: 0.0497
summary(model4) # sig 3.39e-08 -- 30 sep 2021: 1 e -14
summary(model5) # sig 0.0481-- 30 sep 2021:  0.664

summary(model1s) # 0.03 mod, 1 e -10 spp
summary(model2s) # 0.0015 mod, e-16 spp
summary(model3s) # 0.0275 mod, e -07 spp
summary(model4s) # mod and spp both e-16
summary(model5s) # 0.674 mod 0.997 spp


TukeyHSD(model1) # none, lowest is 0.1452022-- 30 sep 2021: none, lowest is 0.13
TukeyHSD(model2) # none, lowest is 0.8256835-- 30 sep 2021: H-B
TukeyHSD(model3) # none, lowest 0.5739302-- 30 sep 2021: none, lowest is 0.06 A-D
TukeyHSD(model4) # mix-iba 0, mix-ibd 0.0002208, mix-ibe 0.0000002, mix-ibh 0.0000002-- 30 sep 2021: H-A, M-A, H-B, M-B, H-D, M-D, H-E, M-E, M-H
TukeyHSD(model5) # none lowest is ibh-ibd 0.0767787-- 30 sep 2021: none, lowest is 0.58

TukeyHSD(model1s) # mod lowest 0.08, 
TukeyHSD(model2s) # 
TukeyHSD(model3s) # 
TukeyHSD(model4s) #  
TukeyHSD(model5s) # 

## aov model -- TRIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL3!=""]) ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])

model1s=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""] +  df$SPECIES[df$MOSTA.1MODEL3!=""])
model2s=aov(df$MEAN_FST[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""]+  df$SPECIES[df$MOSTA.1MODEL3!=""])
model3s=aov(df$MEAN_DXY[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""]+  df$SPECIES[df$MOSTA.1MODEL3!=""])
model4s=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""]+  df$SPECIES[df$MOSTA.1MODEL3!=""])
model5s=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL3!=""]) ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""]+  df$SPECIES[df$MOSTA.1MODEL3!=""])


summary(model1) # sig 0.0472
summary(model2) # not 0.499
summary(model3) # nearly 0.0609
summary(model4) # sig 3.3e-12
summary(model5) # not sig 0.315

summary(model1s) # sig 0.0472
summary(model2s) # not 0.499
summary(model3s) # nearly 0.0609
summary(model4s) # sig 3.3e-12
summary(model5s) # not sig 0.315

TukeyHSD(model1) # sig ibe-iba 0.0402975
TukeyHSD(model2) # none
TukeyHSD(model3) # nearly ibe-iba 0.0676740
TukeyHSD(model4) # sig ibh-iba 0.0000728, mix-iba 0.0000000, ibh-ibe 0.0017143, 
# mix-ibe 0.0000000, mix-ibh 0.0024778
TukeyHSD(model5) # none

TukeyHSD(model1s) # sig ibe-iba 0.0402975
TukeyHSD(model2s) # none
TukeyHSD(model3s) # nearly ibe-iba 0.0676740
TukeyHSD(model4s) # sig ibh-iba 0.0000728, mix-iba 0.0000000, ibh-ibe 0.0017143, 
# mix-ibe 0.0000000, mix-ibh 0.0024778
TukeyHSD(model5s) # none


png("four_panel_figure_chapter_3_nothypotheses.png",height=4,width=6,units = "in",
    res=300)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,7,3,4)],"grey",brewer.pal(8,"Dark2")[c(7)])
## four panels to see if same
par(mfrow=c(2,2))
par(mar=c(4,4,0,0))
boxplot(df$MEAN_FST[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="Mean Fst (cut off)",ylim=c(0,0.21),
        xlab="Best Model",las=2)
#legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
#       col=cols,fill=cols)
boxplot(df$MEAN_DXY[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="Mean Dxy",
        xlab="Best Model",las=2)
boxplot(df$PERCENT_MISSING[df$MOST.1MODEL2!=""]*100~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="Percent Missing",
        xlab="Best Model",las=2)
boxplot(df$MEAN_RECOMB[df$MOST.1MODEL2!=""]/(1e-9)~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="Mean Recombination Rate (x 1^-9)",
        xlab="Best Model",las=2)
dev.off()

## aov model -- bivariate no hypotheses
model1=aov(df$MEAN_RECOMB[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
model2=aov(df$MEAN_FST[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
MODEL3=aov(df$MEAN_DXY[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
model4=aov(df$PERCENT_MISSING[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])

summary(model1) # almost
summary(model2) # not
summary(MODEL3) # sig
summary(model4) # sig

TukeyHSD(model1) # none
TukeyHSD(model2) # none
TukeyHSD(MODEL3) # ibd-abun, ibd-env, ibd-lgm, ibd-mix
TukeyHSD(model4) # mix-abun, mix-env, mix-ibd, mix-lgm, mix-pres, lgm-env

png("missing_boxplot_model2.png")
boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"]*100
        ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"],
        col=cols,
        ylab="Percent Missing",
        xlab="Best Model",las=2)
dev.off()

model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"] 
           ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"])
summary(model4) ## still sig
TukeyHSD(model4) # 


## do a pca?
df_numb = df[,c("PERCENT_MISSING","MEAN_FST","MEAN_DXY","MEAN_RECOMB")]
df_meta = df[,c("DATASET","SPECIES","MOST.1MODEL2","MOSTA.1MODEL2")]

df_meta = df_meta[complete.cases(df_numb),]
df_numb = df_numb[complete.cases(df_numb),]

pca = prcomp(df_numb,center = T,scale. = T)
summary(pca)
df_pca = cbind(df_meta,pca$x)
palette(c(brewer.pal(8,"Dark2")[c(1,7,3,4)],"grey"))
png("pca_plot_by_model.png")
plot(df_pca$PC2,df_pca$PC1,col=as.numeric(as.factor(df_pca$MOSTA.1MODEL2)),
     pch=as.numeric(as.factor(df_pca$MOSTA.1MODEL2)),ylim=c(-3,2))
dev.off()

## look at tajima's 
png("tajimas_boxplots.png")
par(mfrow=c(1,3),mar=c(4,4,0,0))
cols=c(brewer.pal(8,"Dark2")[c(1,7,3,4)],"grey",brewer.pal(8,"Dark2")[c(7)])
boxplot(df$TAJIMAS[df$MOST.1MODEL1!=""]~df$MOST.1MODEL1[df$MOST.1MODEL1!=""],
        col=cols,
        ylab="Tajima's D",
        xlab="",las=2)
boxplot(df$TAJIMAS[df$MOST.1MODEL2!=""]~df$MOST.1MODEL2[df$MOST.1MODEL2!=""],
        col=cols,
        ylab="",
        xlab="Best Model",las=2)
boxplot(df$TAJIMAS[df$MOST.1MODEL3!=""]~df$MOST.1MODEL3[df$MOST.1MODEL3!=""],
        col=cols,
        ylab="",
        xlab="",las=2)
dev.off()

png("tajimas_boxplots_lines.png")
par(mfrow=c(1,3),mar=c(4,4,0,0))
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
boxplot(df$TAJIMAS[df$MOSTA.1MODEL1!=""]~df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""],
        col=cols,
        ylab="Tajima's D",
        xlab="",las=2)
abline(h=0)
boxplot(df$TAJIMAS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,
        ylab="",
        xlab="Best Model",las=2)
abline(h=0)
boxplot(df$TAJIMAS[df$MOSTA.1MODEL3!=""]~df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""],
        col=cols,
        ylab="",
        xlab="",las=2)
abline(h=0)
dev.off()

modelA=aov(df$TAJIMAS[df$MOST.1MODEL1!=""] ~ df$MOST.1MODEL1[df$MOST.1MODEL1!=""])
modelB=aov(df$TAJIMAS[df$MOST.1MODEL2!=""] ~ df$MOST.1MODEL2[df$MOST.1MODEL2!=""])
modelC=aov(df$TAJIMAS[df$MOST.1MODEL3!=""] ~ df$MOST.1MODEL3[df$MOST.1MODEL3!=""])

summary(modelA) # sig
summary(modelB) # sig 0.0107
summary(modelC) # sig

TukeyHSD(modelA) # ibd-abun
TukeyHSD(modelB) # ibd-abun
TukeyHSD(modelC) # env-abun nearly, mix-abun

## negative tajd = lots of rare alleles, recent sweep, expansion after bottleneck, or linkage
## positive tajd = few rares alleles, balancing selection, sudden contraction


cor((df[,c("TAJIMAS","PERCENT_MISSING","MEAN_FST","MEAN_DXY","MEAN_RECOMB")]),
    use = "pairwise.complete.obs")




## chrom lengths
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
png("chrom_len_vs_model_log_notlot.png")
par(mfrow=c(1,2))
boxplot(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],col=cols)
boxplot(log10(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""])~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],col=cols)
#plotrix::gap.boxplot(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],col=cols,
#                     gap=list(top=c(16000,111000),bottom=c(NA,NA)))
title(ylab="Num Scaffolds")
dev.off()

#boxplot(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]*100000)
#        ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
#        col=cols,xlab="Best Model",ylab="Log Chromosome Length")

df_nogen = df[df$DATASET!="GENOME",]

modelA=aov(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelA) # not sig
TukeyHSD(modelA) # 

modelB=aov(df_nogen$CHROM_LENGTH[df_nogen$MOSTA.1MODEL2!=""]
           ~df_nogen$MOSTA.1MODEL2[df_nogen$MOSTA.1MODEL2!=""])
summary(modelB) # not sig
TukeyHSD(modelB) # 

modelC=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]*100000)
           ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
summary(modelC) # not sig
TukeyHSD(modelC) # 


## looking at new stats


pdf("gdm_bivariate_all_sumstats.pdf",height=4,width=6)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
## four panels to see if same
par(mfrow=c(2,2),cex.axis=1)
par(mar=c(0.1,4,0.3,0.1))
#boxplot(df$MEAN_FST[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
#        col=cols,xaxt="n",
#        ylab="Mean Fst",
#        xlab="Best Model",las=2)
plotrix::gap.boxplot(df$MEAN_FST[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
                     gap=list(top=c(0.22,0.56),bottom=c(NA,NA)),
                     axes=F,ylim=c(0.01,0.58),ylab="Mean Fst",xlab="Best Model",
                     xaxt="n",col=cols)
legend("topright",legend=c("IBA","IBD","IBE","IBH","MIX"),
       col=cols,fill=cols)
title(ylab="Mean Fst")
axis(2,labels=c(seq(0,0.2,0.05),0.57),
     at=c(seq(0,0.2,0.05),0.235),tick=T,las=1)
boxplot(df$MEAN_DXY[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Mean Dxy",
        xlab="Best Model",las=1)
boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""]*100~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Percent Missing",
        xlab="Best Model",las=1)
boxplot(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""]/(1e-9)~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Mean Recombination Rate",
        xlab="Best Model",las=2)
boxplot(df$number.fst.peaks[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="N FST Peaks",
        xlab="Best Model",las=2)
boxplot(df$number.fst.peaks[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,500),
        ylab="N FST Peaks (Cropped)",
        xlab="Best Model",las=2)
boxplot(df$number.dxy.peaks[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="N DXY Peaks",
        xlab="Best Model",las=2)
boxplot(df$number.dxy.peaks[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,500),
        ylab="N DXY Peaks (Cropped)",
        xlab="Best Model",las=2)
boxplot(df$Number.islands[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="N Islands",
        xlab="Best Model",las=2)
boxplot(df$TAJIMAS[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Tajima's D",
        xlab="Best Model",las=2)
boxplot(df$number.sweeps[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="N Sweeps",
        xlab="Best Model",las=2)
boxplot(df$number.sweeps[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,100),
        ylab="N Sweeps (Cropped)",
        xlab="Best Model",las=2)
boxplot(df$number.dxy.lows[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="N Dxy Lows",
        xlab="Best Model",las=2)
boxplot(df$number.dxy.lows[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,500),
        ylab="N Dxy Lows (Cropped)",
        xlab="Best Model",las=2)
boxplot(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,20000),
        ylab="Chrom Length (not Genome)",
        xlab="Best Model",las=2)
boxplot(df$prop.fst.peak[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Proportion FST Peak",
        xlab="Best Model",las=2)
boxplot(df$prop.dxy.peak[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Proportion DXY Peak",
        xlab="Best Model",las=2)
boxplot(df$prop.isl[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Proportion Island",
        xlab="Best Model",las=2)
boxplot(df$prop.swp[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Proportion Sweep",
        xlab="Best Model",las=2)
boxplot(df$prop.swp[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,0.01),
        ylab="Proportion Sweep (Cropped)",
        xlab="Best Model",las=2)
boxplot(df$prop.dxy.low[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",
        ylab="Proportion Low DXY",
        xlab="Best Model",las=2)
boxplot(df$prop.dxy.low[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xaxt="n",ylim=c(0,0.35),
        ylab="Proportion Low DXY (Cropped)",
        xlab="Best Model",las=2)

dev.off()


## THIS SECTION NOT WORKING BECAUSE USES OLD NOTATION
small=df[df$MOSTA.1MODEL2!="",c("MOSTA.1MODEL2","TAJIMAS","proportion_fst_outliers_under_5sd","proportion_fst_outliers_over_5sd","prop.swp","prop.isl","prop.fst.peak","prop.dxy.peak","prop.dxy.low","PERCENT_MISSING","number.sweeps","Number.islands","number.fst.peaks","number.dxy.peaks","number.dxy.lows","number_fst_outliers_under_5sd","number_fst_outliers_over_5sd","MEAN_RECOMB","MEAN_FST","MEAN_DXY","HZAR_width","HZAR_center","CHROM_LENGTH")]
agg=aggregate(cbind(small[,2],small[,3],small[,4],small[,5],small[,6],
                    small[,7],small[,8],small[,9],
                    small[,10],small[,11],small[,12],
                    small[,13],small[,14],small[,15],small[,16],
                    small[,17],small[,18],small[,19],small[,20],
                    small[,21],small[,22],small[,23])
              ~small[,1],FUN=function(x){mean(x,na.rm=T)})
colnames(agg) = colnames(small)#[c(17,1:16)]

pdf("gdm_bivariate_all_barplots.pdf",height=4,width=6)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
## four panels to see if same
par(mfrow=c(2,2),cex.axis=1)
par(mar=c(0.1,4,0.3,0.1))
for(i in 2:17){
  print(i)
  barplot(agg[,i],col=cols,ylab=colnames(agg)[i])
}
dev.off()

