df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate/bivariate_gdm_results.csv",
                sep="\t",
                header=T,
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
plot(numonly_nofstout)

cor(numonly,use="pairwise.complete.obs")
cor(numonly_nofstout,use="pairwise.complete.obs")

par(mfrow=c(1,2))
corrplot::corrplot(cor(numonly,use="pairwise.complete.obs"),
                   method = "number",diag=T,type="upper")
corrplot::corrplot(cor(numonly_nofstout,use="pairwise.complete.obs"),
                   method="number",diag=T,type="upper")

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






plot(df$MEAN_FST,df$PERCENT_MISSING)
plot(df$MEAN_DXY,df$PERCENT_MISSING)
plot(df$MEAN_RECOMB,df$PERCENT_MISSING)

cor(df$MEAN_FST,df$PERCENT_MISSING,use="pairwise.complete.obs") # -0.337
cor(df$MEAN_DXY,df$PERCENT_MISSING,use="pairwise.complete.obs") # 0.267
cor(df$MEAN_RECOMB,df$PERCENT_MISSING,use="pairwise.complete.obs") # 0.08

library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
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

model=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
summary(model) 
## sig uni, mixed different
## not significant biv
## almost sig multim mixed-inh almost sig
TukeyHSD(model)

agg=aggregate(df$PERCENT_MISSING~df$SPECIES,FUN=function(x){sd(x,na.rm=T)})


png("four_panel_figure_chapter_3_ibd2mix.png",height=4,width=6,units = "in",
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
dev.off()
## sig uni, mixed different
## not significant biv
## almost sig multim mixed-inh almost sig


png("four_panel_figure_chapter_3_univariate_ibd2mix.png",height=4,width=6,units = "in",
    res=300)
library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,3,7,4)],"grey")
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



df = df[df$SPECIES!="NITENS",]

## aov model -- BIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL2!=""] ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]) ~ df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""])

summary(model1) # sig 0.0498
summary(model2) # not
summary(model3) # sig 0.0245
summary(model4) # sig 1.96e-9 -- still sig without nitens
summary(model5) # not sig 0.115

TukeyHSD(model1) # nearly ibe-iba 0.0565688
TukeyHSD(model2) # none
TukeyHSD(model3) # sig ibd-iba 0.0340591
TukeyHSD(model4) # ibd-iba 0.0045044 mix-iba 0 mix-ibd 0.0008079 ibh-ibe 0.0042963 mix-ibe 0 mix-ibh 0.0127600
## ibh comparisons not sig without nitens 
TukeyHSD(model5) # ibd-iba 0.036m ibe-ibd 0.066, ibh-ibd = 0.027



## aov model -- UNIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL1!=""] ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL1!=""]) ~ df$MOSTA.1MODEL1[df$MOSTA.1MODEL1!=""])


summary(model1) # sig 0.0455
summary(model2) # not 0.846
summary(model3) # not 0.527
summary(model4) # sig 3.39e-08
summary(model5) # sig 0.0481

TukeyHSD(model1) # none, lowest is 0.1452022
TukeyHSD(model2) # none, lowest is 0.8256835
TukeyHSD(model3) # none, lowest 0.5739302
TukeyHSD(model4) # mix-iba 0, mix-ibd 0.0002208, mix-ibe 0.0000002, mix-ibh 0.0000002
TukeyHSD(model5) # none lowest is ibh-ibd 0.0767787


## aov model -- TRIVARIATE
model1=aov(df$MEAN_RECOMB[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model2=aov(df$MEAN_FST[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model3=aov(df$MEAN_DXY[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model4=aov(df$PERCENT_MISSING[df$MOSTA.1MODEL3!=""] ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])
model5=aov(log(df$CHROM_LENGTH[df$MOSTA.1MODEL3!=""]) ~ df$MOSTA.1MODEL3[df$MOSTA.1MODEL3!=""])

summary(model1) # sig 0.0472
summary(model2) # not 0.499
summary(model3) # nearly 0.0609
summary(model4) # sig 3.3e-12
summary(model5) # not sig 0.315

TukeyHSD(model1) # sig ibe-iba 0.0402975
TukeyHSD(model2) # none
TukeyHSD(model3) # nearly ibe-iba 0.0676740
TukeyHSD(model4) # sig ibh-iba 0.0000728, mix-iba 0.0000000, ibh-ibe 0.0017143, 
# mix-ibe 0.0000000, mix-ibh 0.0024778
TukeyHSD(model5) # none



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

boxplot(df$PERCENT_MISSING[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"]*100
        ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!="" & df$SPECIES!="NITENS"],
        col=cols,
        ylab="Percent Missing",
        xlab="Best Model",las=2)

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
plot(df_pca$PC2,df_pca$PC1,col=as.numeric(as.factor(df_pca$MOSTA.1MODEL2)),
     pch=as.numeric(as.factor(df_pca$MOSTA.1MODEL2)),ylim=c(-3,2))


## look at tajima's 
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
plotrix::gap.boxplot(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],col=cols,
                     gap=list(top=c(16000,111000),bottom=c(NA,NA)))
title(ylab="Num Scaffolds")

boxplot(log(df$CHROM_LENGTH[df$MOSTA.1MODEL2!=""]*100000)
        ~df$MOSTA.1MODEL2[df$MOSTA.1MODEL2!=""],
        col=cols,xlab="Best Model",ylab="Log Chromosome Length")

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
