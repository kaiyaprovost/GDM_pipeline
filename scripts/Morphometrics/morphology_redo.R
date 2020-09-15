file = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"

agg = read.table(file,sep="\t",header=T,stringsAsFactors = F,fill=T)
agg = unique(agg)

head(agg)
numagg = agg[,c("BILL.HEIGHT","BILL.LENGTH","BILL.WIDTH","TARSUS.LENGTH",
                "WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES",
                "TAIL","KIPPSINDEX_CALC","BEAKBASEAREA_CALC",
                "BEAKVOL_CALC","BEAKLATERALSURFACE_CALC","BEAKTOTALSURFACE_CALC",
                "BEAKAREAVOLUME_CALC","BEAKTOTAREAVOLUME_CALC")]

agg = agg[complete.cases(numagg),c("BILL.HEIGHT","BILL.LENGTH","BILL.WIDTH","TARSUS.LENGTH",
                                   "WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES",
                                   "TAIL","KIPPSINDEX_CALC","BEAKBASEAREA_CALC",
                                   "BEAKVOL_CALC","BEAKLATERALSURFACE_CALC","BEAKTOTALSURFACE_CALC",
                                   "BEAKAREAVOLUME_CALC","BEAKTOTAREAVOLUME_CALC",
                                   "LAT","LONG","GENETIC.SIDE","SPP","WHICH.SIDE.OF.CFB")]
numagg = numagg[complete.cases(numagg),]
agg=unique(agg)

pca = prcomp(agg[,c("BILL.HEIGHT","BILL.LENGTH","BILL.WIDTH","TARSUS.LENGTH",
                    "WING.LENGTH.TO.PRIMARIES","WING.LENGTH.TO.SECONDARIES",
                    "TAIL","KIPPSINDEX_CALC","BEAKBASEAREA_CALC",
                    "BEAKVOL_CALC","BEAKLATERALSURFACE_CALC","BEAKTOTALSURFACE_CALC",
                    "BEAKAREAVOLUME_CALC","BEAKTOTAREAVOLUME_CALC")],
             center = T,scale. = T)
pca$rotation[,1:3]
summary(pca)

dat = as.data.frame(pca$x)
axes_full = cbind(dat,agg)
plot(axes_full$PC1,dat$PC2,col=as.numeric(as.factor(axes_full$SPP)))
plot(axes_full$PC1,dat$PC3,col=as.numeric(as.factor(axes_full$SPP)))

boxplot((axes_full$PC1/max(abs(axes_full$PC1)))~axes_full$SPP,
        ylab="relative PC value",ylim=c(-1,1),at=seq(1,30,3),
        xaxt="n",xlab="",col="white",xlim=c(1,30))
boxplot((axes_full$PC2/max(abs(axes_full$PC2)))~axes_full$SPP,
        ylab="relative PC value",at=seq(2,30,3),add=T,las=2,
        yaxt="n",col="red")
boxplot((axes_full$PC3/max(abs(axes_full$PC3)))~axes_full$SPP,
        ylab="relative PC value",at=seq(3,30,3),add=T,xaxt="n",
        yaxt="n",col="cyan")
abline(v=seq(0.5,30.5,3),lty=2,col="grey")

png("all_pcs_plot.png",width=900,height=400); {
  par(mfrow=c(1,4),pt.cex=2,cex=2
      #,mar=c(2,2,0,0)
  )
  m <- rbind(c(1,1,1), c(2,3,4), c(2,3,4), c(2,3,4), c(2,3,4), c(2,3,4))
  #print(m)
  layout(m)
  #layout.show(4)
  palette(c("black","orange","goldenrod","green","pink",
            "magenta","purple","blue","grey","red"))
  par(mar=c(0,0,0,0))
  plot.new()
  legend("center",legend=sort(unique(axes_full$SPP)),ncol=5,bty="n",cex=2,pt.cex=2,
         col=1:10,pch=1:10)
  #box()
  par(mar=c(4,4,0,0))
  plot(axes_full$PC1,axes_full$PC2,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),ylab="PC2",xlab="PC1",type="n")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC1","PC2")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC1, toplot$PC2, levels=c(0.75),lty=1,pch=i,
                     add=T,col=i,plot.points=T,center.pch = F)
  }
  
  plot(axes_full$PC3,axes_full$PC2,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),ylab="PC2",xlab="PC3",type="n")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC3","PC2")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC3, toplot$PC2, levels=c(0.75),lty=1,pch=i,
                     add=T,col=i,plot.points=T,center.pch = F)
  }
  
  
  plot(axes_full$PC1,axes_full$PC3,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),ylab="PC3",xlab="PC1",type="n")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC1","PC3")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC1, toplot$PC3, levels=c(0.75),lty=1,pch=i,
                     add=T,col=i,plot.points=T,center.pch = F)
  }
  
  
}; dev.off()


pdf("pca_presentation.pdf",height=4,width=7); {
  par(mfrow=c(1,2),bg="black",col.axis="white",
      col.lab="white",col="white")
  palette(c("white","orange","goldenrod","green","pink",
            "magenta","purple","cyan","grey","red"))
  par(mar=c(2,2,0,0))
  plot(axes_full$PC1,axes_full$PC2,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),ylab="PC2",xlab="PC1",type="n")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC1","PC2")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC1, toplot$PC2, levels=c(0.75),lty=1,
                     add=T,col=i,plot.points=T,center.pch = F)
  }
  plot(axes_full$PC3,axes_full$PC2,col=as.factor(axes_full$SPP),
       pch=as.numeric(as.factor(axes_full$SPP)),ylab="PC2",xlab="PC3",type="n")
  for(i in 1:length(unique(axes_full$SPP))) {
    spp = sort(unique(axes_full$SPP))[i]
    toplot = axes_full[axes_full$SPP==spp,c("PC3","PC2")]
    toplot = toplot[complete.cases(toplot),]
    car::dataEllipse(toplot$PC3, toplot$PC2, levels=c(0.75),lty=1,
                     add=T,col=i,plot.points=T,center.pch = F)
  }
}; dev.off()




library(dabestr)
print("PC1"); {
  temp = axes_full
  temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
  temp = temp[,c("GROUP","PC1")]
  temp = temp[complete.cases(temp), ]
  bav = dabest(temp,GROUP,PC1, idx=(list(c("BIL_SON","BIL_CHI"),
                                          c("BRU_SON","BRU_CHI"),
                                          c("MEL_SON","MEL_CHI"),
                                          c("NIT_SON","NIT_CHI"),
                                          c("SIN_SON","SIN_CHI"),
                                          c("BEL_SON","BEL_CHI"),
                                          c("CRI_SON","CRI_CHI"),
                                          c("CUR_SON","CUR_CHI"),
                                          c("FLA_SON","FLA_CHI"),
                                          c("FUS_SON","FUS_CHI"))))
  g=plot(x=bav,group.summaries=NULL,
         rawplot.ylabel="PC1",color.column="GROUP",
         #palette="scale_colour_grey",
         #theme=theme_transparentwhite(),
         effsize.markersize=2)
  print(g)
  ggplot2::ggsave("DABEST_PC1.png",
         g,bg="transparent",width=23,height=7,units="in")
}
print("PC2"); {
  temp = axes_full
  temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
  temp = temp[,c("GROUP","PC2")]
  temp = temp[complete.cases(temp), ]
  bav = dabest(temp,GROUP,PC2, idx=(list(c("BIL_SON","BIL_CHI"),
                                         c("BRU_SON","BRU_CHI"),
                                         c("MEL_SON","MEL_CHI"),
                                         c("NIT_SON","NIT_CHI"),
                                         c("SIN_SON","SIN_CHI"),
                                         c("BEL_SON","BEL_CHI"),
                                         c("CRI_SON","CRI_CHI"),
                                         c("CUR_SON","CUR_CHI"),
                                         c("FLA_SON","FLA_CHI"),
                                         c("FUS_SON","FUS_CHI"))))
  g=plot(x=bav,group.summaries=NULL,
         rawplot.ylabel="PC2",color.column="GROUP",
         #palette="scale_colour_grey",
         #theme=theme_transparentwhite(),
         effsize.markersize=2)
  print(g)
  ggplot2::ggsave("DABEST_PC2.png",
                  g,bg="transparent",width=23,height=7,units="in")
}
print("PC3"); {
  temp = axes_full
  temp$GROUP = paste(substr(temp$SPP,1,3),substr(temp$WHICH.SIDE.OF.CFB,1,3),sep="_")
  temp = temp[,c("GROUP","PC3")]
  temp = temp[complete.cases(temp), ]
  bav = dabest(temp,GROUP,PC3, idx=(list(c("BIL_SON","BIL_CHI"),
                                         c("BRU_SON","BRU_CHI"),
                                         c("MEL_SON","MEL_CHI"),
                                         c("NIT_SON","NIT_CHI"),
                                         c("SIN_SON","SIN_CHI"),
                                         c("BEL_SON","BEL_CHI"),
                                         c("CRI_SON","CRI_CHI"),
                                         c("CUR_SON","CUR_CHI"),
                                         c("FLA_SON","FLA_CHI"),
                                         c("FUS_SON","FUS_CHI"))))
  g=plot(x=bav,group.summaries=NULL,
         rawplot.ylabel="PC3",color.column="GROUP",
         #palette="scale_colour_grey",
         #theme=theme_transparentwhite(),
         effsize.markersize=2)
  print(g)
  ggplot2::ggsave("DABEST_PC3.png",
                  g,bg="transparent",width=23,height=7,units="in")
}

