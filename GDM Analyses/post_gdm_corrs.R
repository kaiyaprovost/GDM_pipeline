# post gdm assessment of env PCs against the morphologies 

library(gtools)
library(raster)

morph = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
morphdf = read.csv(morph)
morphdf=unique(morphdf)
morphdf$CATALOG.NUMBER = stringr::str_replace_all(morphdf$CATALOG.NUMBER,"[^[:alnum:]]","")
morphdf = morphdf[morphdf$CATALOG.NUMBER!="",]

morphpca = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019_PCASCORES_FULL.csv"
pcadf = read.csv(morphpca)
pcadf=unique(pcadf)
pcadf$CATALOG.NUMBER = stringr::str_replace_all(pcadf$CATALOG.NUMBER,"[^[:alnum:]]","")

alldf = merge(morphdf,pcadf,all=T)

preclayers = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/Deserts_Bioclim_PrecPCA_From_GDM.tif")
templayers = stack("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/Deserts_Bioclim_TempPCA_From_GDM.tif")

alllocs = alldf[,c("LAT","LONG")]


precdata = as.data.frame(extract(preclayers,alllocs[,2:1]))
tempdata = as.data.frame(extract(templayers,alllocs[,2:1]))

colnames(precdata) = c("PC1P","PC2P","PC3P","PC4P","PC5P")
colnames(tempdata) = c("PC1T","PC2T","PC3T","PC4T","PC5T")

envdata = cbind(precdata,tempdata)
envdata$CATALOG.NUMBER = alldf$CATALOG.NUMBER

locsenv = merge(alldf,envdata,all=T)

correlation=cor(locsenv[,c(2:10,14:20,48:71)],use="pairwise.complete.obs")
smallcorr = t(correlation[c(1:7,10:19),c(8:9,31:33,36:38)])

png("AllSpp_corrplot_env_morph.png")
corrplot::corrplot(round(smallcorr,1),is.corr=F,method="number",main="\nAll Spp",cl.lim=c(-1,1))
dev.off()

for (spp in unique(locsenv$SPP)) {
  print(spp)
  if(spp != "CARDINALIS") {
    locssmall = locsenv[locsenv$SPP==spp,]
    nrow(locssmall)
    correlation=cor(locssmall[,c(2:10,14:20,48:71)],use="pairwise.complete.obs")
    smallcorr = t(correlation[c(1:7,10:19),c(8:9,31:33,36:38)])
    
    png(paste(spp,"_corrplot_env_morph.png",sep=""))
    corrplot::corrplot(round(smallcorr,1),is.corr=F,method="number",main=paste("\n",spp,"N=",nrow(locssmall)),cl.lim=c(-1,1))
    dev.off()
    
  }
}

## real quick correlation bwtween the pctemp and pcprec
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Precipitation_Temperature_PCA_correlation.png")
corrplot::corrplot(cor(locsenv[,c(62:64,67:69)],
                       use = "pairwise.complete.obs"),
                   method="color")
dev.off()

## STRUCTURE
## do boxplot of which side of barrier to morph pca
sppnames =c("bellii",
            "bilineata",
            "brunneicapillus",
            "crissale",
            "curvirostre",
            "flaviceps",
            "fusca",
            "melanura",
            "nitens",
            "sinuatus")

structure = unique(pcadf[,c(1:3,29:34)])
structure = structure[complete.cases(structure),]
structure$SPP = substr(structure$SPP,1,3)
structure$WHICH.SIDE.OF.CFB = substr(structure$WHICH.SIDE.OF.CFB,1,3)
structure$WHICH.SIDE.OF.CFB[structure$WHICH.SIDE.OF.CFB=="UNC"] = "HYB" 

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_MorphPC1_Structure2.png",
    height=450,width=1000); {
      par(mfrow=c(2,5),mar=c(4,4,2,1))
      for (species in unique(sort(structure$SPP))) {
        print(species)
        numcats = length(unique(structure$WHICH.SIDE.OF.CFB[structure$SPP==species]))
        if(numcats > 2) {
          colors=c("cyan","yellow","green")
          names = c("Chi.","Hybrid?","Son.")
        } else {
          colors = c("cyan","green")
          names = c("Chi.","Son.")
        }
        
        boxplot(structure$PC1[structure$SPP==species] ~ structure$WHICH.SIDE.OF.CFB[structure$SPP==species],
                las=1,horizontal=F,main=species,
                xlab="Side of Barrier",ylab="PC1 Morphology",cex=1,col=colors,
                names=names
        )
      }
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_MorphPC2_Structure2.png",
    height=450,width=1000); {
      par(mfrow=c(2,5),mar=c(4,4,2,1))
      for (species in unique(sort(structure$SPP))) {
        print(species)
        numcats = length(unique(structure$WHICH.SIDE.OF.CFB[structure$SPP==species]))
        if(numcats > 2) {
          colors=c("cyan","yellow","green")
          names = c("Chi.","Hybrid?","Son.")
        } else {
          colors = c("cyan","green")
          names = c("Chi.","Son.")
        }
        
        boxplot(structure$PC2[structure$SPP==species] ~ structure$WHICH.SIDE.OF.CFB[structure$SPP==species],
                las=1,horizontal=F,main=species,
                xlab="Side of Barrier",ylab="PC1 Morphology",cex=1,col=colors,
                names=names
        )
      }
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_MorphPC3_Structure2.png",
    height=450,width=1000); {
      par(mfrow=c(2,5),mar=c(4,4,2,1))
      for (species in unique(sort(structure$SPP))) {
        print(species)
        numcats = length(unique(structure$WHICH.SIDE.OF.CFB[structure$SPP==species]))
        if(numcats > 2) {
          colors=c("cyan","yellow","green")
          names = c("Chi.","Hybrid?","Son.")
        } else {
          colors = c("cyan","green")
          names = c("Chi.","Son.")
        }
        
        boxplot(structure$PC3[structure$SPP==species] ~ structure$WHICH.SIDE.OF.CFB[structure$SPP==species],
                las=1,horizontal=F,main=species,
                xlab="Side of Barrier",ylab="PC1 Morphology",cex=1,col=colors,
                names=names
        )
      }
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_MorphPC_Structure.png",
    height=480,width=480*3); {
      par(mfrow=c(1,3))
      
      #png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_str.png")
      boxplot(structure$PC1 ~ structure$WHICH.SIDE.OF.CFB + structure$SPP,las=2,horizontal=T,
              xlab="PC1 (Size)",cex=1,col=c("cyan","yellow","green"),ylab="",yaxt="n")
      axis(2,at=c(seq(2,29,3)),labels=substr(sppnames,1,3),las=2)
      lapply(seq(0.5,33.5,3),FUN=function(x) {abline(h=x,col="grey",lty=3)})
      legend("bottomleft",legend=c("Sonoran","Intermediate","Chihuahuan"),
             fill=c("green","yellow","cyan"))
      #dev.off()
      
      #png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_str.png")
      boxplot(structure$PC2 ~ structure$WHICH.SIDE.OF.CFB + structure$SPP,las=2,horizontal=T,
              xlab="PC2 (Beak Shape)",cex=1,col=c("cyan","yellow","green"),ylab="",yaxt="n")
      axis(2,at=c(seq(2,29,3)),labels=substr(sppnames,1,3),las=2)
      lapply(seq(0.5,33.5,3),FUN=function(x) {abline(h=x,col="grey",lty=3)})
      legend("bottomright",legend=c("Sonoran","Intermediate","Chihuahuan"),
             fill=c("green","yellow","cyan"))
      #dev.off()
      
      #png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_str.png")
      boxplot(structure$PC3 ~ structure$WHICH.SIDE.OF.CFB + structure$SPP,las=2,horizontal=T,
              xlab="PC3 (Wing Shape)",cex=1,col=c("cyan","yellow","green"),ylab="",yaxt="n",
              ylim=c(-4,4))
      axis(2,at=c(seq(2,29,3)),labels=substr(sppnames,1,3),las=2)
      lapply(seq(0.5,33.5,3),FUN=function(x) {abline(h=x,col="grey",lty=3)})
      legend("topleft",legend=c("Sonoran","Intermediate","Chihuahuan"),
             fill=c("green","yellow","cyan"))
    }
dev.off()

## environment
environment = locsenv[,c(9:10,13,48:50,62:64,67:69)]
environment = environment[complete.cases(environment),]

## pc1
collist = c("red","goldenrod","orange","green","blue","cyan","magenta","purple","black","brown","grey")
palette(collist)
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_Precip.png",
    height=480,width=480*3); {
      
      par(mfrow=c(1,3))
      
      plot(environment$PC1P,environment$PC1,xlab="PC1Precipitation",ylab="Morphology PC1",
           col=as.numeric(as.factor(as.character(environment$SPP))),
           pch=as.numeric(as.factor(as.character(environment$SPP))))
      for(species in sort(unique(as.factor(as.character(environment$SPP))))) {
        num = which(species == sort(unique(as.factor(as.character(environment$SPP)))))
        color = collist[num]
        model = lm(environment$PC1[environment$SPP==species]~environment$PC1P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,
                 col=color,lty=num%%3+1,
                 lwd=3)
        } else {
          abline(model,
                 col=color,lty=num%%3+1,
                 lwd=1)
        }
      }
      legend("topleft",legend=substr(sort(unique(as.factor(as.character(environment$SPP)))),1,3),
             col=collist[1:10],pch=1:10,cex=1,ncol=2)
      
      
      plot(environment$PC2P,environment$PC1,xlab="PC2Precipitation",ylab="Morphology PC1",
           col=as.numeric(as.factor(as.character(environment$SPP))),
           pch=as.numeric(as.factor(as.character(environment$SPP))))
      for(species in sort(unique(as.factor(as.character(environment$SPP))))) {
        num = which(species == sort(unique(as.factor(as.character(environment$SPP)))))
        color = collist[num]
        model = lm(environment$PC1[environment$SPP==species]~environment$PC2P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,
                 col=color,lty=num%%3+1,
                 lwd=3)
        } else {
          abline(model,
                 col=color,lty=num%%3+1,
                 lwd=1)
        }
      }
      legend("topleft",legend=substr(sort(unique(as.factor(as.character(environment$SPP)))),1,3),
             col=collist[1:10],pch=1:10,cex=1,ncol=2)
      
      
      plot(environment$PC3P,environment$PC1,xlab="PC3Precipitation",ylab="Morphology PC1",
           col=as.numeric(as.factor(as.character(environment$SPP))),
           pch=as.numeric(as.factor(as.character(environment$SPP))))
      for(species in sort(unique(as.factor(as.character(environment$SPP))))) {
        num = which(species == sort(unique(as.factor(as.character(environment$SPP)))))
        color = collist[num]
        model = lm(environment$PC1[environment$SPP==species]~environment$PC3P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,
                 col=color,lty=num%%3+1,
                 lwd=3)
        } else {
          abline(model,
                 col=color,lty=num%%3+1,
                 lwd=1)
        }
      }
      legend("topleft",legend=substr(sort(unique(as.factor(as.character(environment$SPP)))),1,3),
             col=collist[1:10],pch=1:10,cex=1,ncol=2)
    }
dev.off()



png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_PC1P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC1P[environment$SPP==species],environment$PC1[environment$SPP==species],
             xlab="PC1Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC1[environment$SPP==species]~environment$PC1P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_PC2P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC2P[environment$SPP==species],environment$PC1[environment$SPP==species],
             xlab="PC2Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC1[environment$SPP==species]~environment$PC2P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_PC3P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC3P[environment$SPP==species],environment$PC1[environment$SPP==species],
             xlab="PC3Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC1[environment$SPP==species]~environment$PC3P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

## now temps
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_PC1T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC1T[environment$SPP==species],environment$PC1[environment$SPP==species],
             xlab="PC1Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC1[environment$SPP==species]~environment$PC1T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_PC2T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC2T[environment$SPP==species],environment$PC1[environment$SPP==species],
             xlab="PC2Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC1[environment$SPP==species]~environment$PC2T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC1M_PC3T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC3T[environment$SPP==species],environment$PC1[environment$SPP==species],
             xlab="PC3Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC1[environment$SPP==species]~environment$PC3T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()



png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_PC1P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC1P[environment$SPP==species],environment$PC2[environment$SPP==species],
             xlab="PC1Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC2[environment$SPP==species]~environment$PC1P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_PC2P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC2P[environment$SPP==species],environment$PC2[environment$SPP==species],
             xlab="PC2Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC2[environment$SPP==species]~environment$PC2P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_PC3P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC3P[environment$SPP==species],environment$PC2[environment$SPP==species],
             xlab="PC3Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC2[environment$SPP==species]~environment$PC3P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

## now temps
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_PC1T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC1T[environment$SPP==species],environment$PC2[environment$SPP==species],
             xlab="PC1Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC2[environment$SPP==species]~environment$PC1T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_PC2T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC2T[environment$SPP==species],environment$PC2[environment$SPP==species],
             xlab="PC2Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC2[environment$SPP==species]~environment$PC2T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC2M_PC3T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC3T[environment$SPP==species],environment$PC2[environment$SPP==species],
             xlab="PC3Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC2[environment$SPP==species]~environment$PC3T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()




png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_PC1P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC1P[environment$SPP==species],environment$PC3[environment$SPP==species],
             xlab="PC1Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC3[environment$SPP==species]~environment$PC1P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_PC2P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC2P[environment$SPP==species],environment$PC3[environment$SPP==species],
             xlab="PC2Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC3[environment$SPP==species]~environment$PC2P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_PC3P.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC3P[environment$SPP==species],environment$PC3[environment$SPP==species],
             xlab="PC3Precipitation",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC3[environment$SPP==species]~environment$PC3P[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

## now temps
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_PC1T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC1T[environment$SPP==species],environment$PC3[environment$SPP==species],
             xlab="PC1Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC3[environment$SPP==species]~environment$PC1T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_PC2T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC2T[environment$SPP==species],environment$PC3[environment$SPP==species],
             xlab="PC2Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC3[environment$SPP==species]~environment$PC2T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/postGDM_relationships_PC3M_PC3T.png",
    height=450,width=1000); {
      
      
      par(mfrow=c(2,5),mar=c(4,4,1,1))
      
      for(species in sort(unique(environment$SPP))) {
        plot(environment$PC3T[environment$SPP==species],environment$PC3[environment$SPP==species],
             xlab="PC3Temperature",ylab="Morphology PC1",
             main=species)
        model = lm(environment$PC3[environment$SPP==species]~environment$PC3T[environment$SPP==species])
        adjr = summary(model)$adj.r.squared
        mtext(paste("AdjRsq:",round(adjr,2)),side=4)
        pval = summary(model)$coefficients[2,4]  
        if(pval < 0.05) {
          abline(model,lwd=3,col="red")
        } else {
          abline(model,lwd=1,lty=3)
        }
      }
      
    }
dev.off()

## ABUN
## this is a resistance layer -- need to figure out how to visualize? 
abunfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",
                       pattern="clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt",
                       recursive = T,
                       full.names = T)

pcfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",
                     pattern="distancematrix_pc",
                     recursive = T,
                     full.names = T)



for (i in 1:length(abunfiles)) {
  print(i)
  file = abunfiles[i]
  pc1file = pcfiles[(((i-1)*3)+1)]
  pc2file = pcfiles[(((i-1)*3)+2)]
  pc3file = pcfiles[(((i-1)*3)+3)]
  species = strsplit(file,split="/")[[1]][11]
  print(species)
  
  png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/",
            species,"_morphology_abundance_distance_plots.png",sep=""),
      height=450,width=1000)
  
  ABUN=read.csv(file,sep="\t")
  if(length(names(ABUN)) == 1) {
    ABUN=read.csv(file)
  }
  PC1 = read.csv(pc1file)
  PC2 = read.csv(pc2file)
  PC3 = read.csv(pc2file)
  
  
  
  colnames(ABUN)[1] = "Sample"
  rownames(ABUN) = stringr::str_replace_all(rownames(ABUN),"[^[:alnum:]]","")
  ABUN$Sample = stringr::str_replace_all(ABUN$Sample,"[^[:alnum:]]","")
  colnames(PC1)[1] = "Sample"
  rownames(PC1) = stringr::str_replace_all(rownames(PC1),"[^[:alnum:]]","")
  PC1$Sample = stringr::str_replace_all(PC1$Sample,"[^[:alnum:]]","")  
  colnames(PC2)[1] = "Sample"
  rownames(PC2) = stringr::str_replace_all(rownames(PC2),"[^[:alnum:]]","")
  PC2$Sample = stringr::str_replace_all(PC2$Sample,"[^[:alnum:]]","")  
  colnames(PC3)[1] = "Sample"
  rownames(PC3) = stringr::str_replace_all(rownames(PC3),"[^[:alnum:]]","")
  PC3$Sample = stringr::str_replace_all(PC3$Sample,"[^[:alnum:]]","")  
  
  colnames(ABUN) = stringr::str_replace_all(colnames(ABUN),"[^[:alnum:]]","")
  colnames(PC1) = stringr::str_replace_all(colnames(PC1),"[^[:alnum:]]","")
  colnames(PC2) = stringr::str_replace_all(colnames(PC2),"[^[:alnum:]]","")
  colnames(PC3) = stringr::str_replace_all(colnames(PC3),"[^[:alnum:]]","")
  
  
  tokeep = intersect(PC1$Sample,ABUN$Sample)
  
  ABUN = ABUN[which(ABUN$Sample %in% tokeep),which(colnames(ABUN) %in% tokeep)]
  PC1 = PC1[which(PC1$Sample %in% tokeep),which(colnames(PC1) %in% tokeep)]
  PC2 = PC2[which(PC2$Sample %in% tokeep),which(colnames(PC2) %in% tokeep)]
  PC3 = PC3[which(PC3$Sample %in% tokeep),which(colnames(PC3) %in% tokeep)]
  
  
  pairwise_abun = data.frame( t(combn(names(ABUN),2)), 
                              dist=t(ABUN)[lower.tri(ABUN)] )
  pairwise_PC1 = data.frame( t(combn(names(PC1),2)),dist=t(PC1)[lower.tri(PC1)] )
  pairwise_PC2 = data.frame( t(combn(names(PC2),2)),dist=t(PC2)[lower.tri(PC2)] )
  pairwise_PC3 = data.frame( t(combn(names(PC3),2)),dist=t(PC3)[lower.tri(PC3)] )
  
  names(pairwise_abun) = c("FIRST","SECOND","ABUNDIST")
  names(pairwise_PC1) = c("FIRST","SECOND","PC1MDIST")
  names(pairwise_PC2) = c("FIRST","SECOND","PC2MDIST")
  names(pairwise_PC3) = c("FIRST","SECOND","PC3MDIST")
  
  mergeddata = merge(pairwise_abun,pairwise_PC1)
  mergeddata = merge(mergeddata,pairwise_PC2)
  mergeddata = merge(mergeddata,pairwise_PC3)
  
  par(mfrow=c(1,3))
  
  for(j in 4:6) {
    plot(mergeddata$ABUNDIST,mergeddata[,j],
         main=species,xlab="Abundance Distance",
         ylab=paste("Morphology PC",j-3," Distance",sep=""))
    model = lm(mergeddata[,j]~mergeddata$ABUNDIST)
    adjr = summary(model)$adj.r.squared
    mtext(paste("AdjRsq:",round(adjr,2)),side=4)
    pval = summary(model)$coefficients[2,4]  
    if(pval < 0.05) {
      abline(model,lwd=3,col="red")
    } else {
      abline(model,lwd=1,lty=3)
    }
  }
  
  dev.off()
}


