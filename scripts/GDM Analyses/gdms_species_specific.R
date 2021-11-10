setwd("~")

df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate_gdm_results_useold.csv",
                sep=",",
                header=T,fill=T,
                stringsAsFactors = F,
                skip = 0)

df=df[df$SPECIES!="",]


library(RColorBrewer)
cols=c(brewer.pal(8,"Dark2")[c(1,2,3,7,4)],"grey")

for(spp in c("SINUATUS")){
#for(spp in unique(df$SPECIES)){
  df2 = df[df$SPECIES==spp,]
  print(spp)
  #pdf(paste(spp,"_boxplots_gdms_each_species.pdf",sep=""))
  par(mfrow=c(1,3))
  # print("FST")
  #   boxplot(df2$MEAN_FST_100[df2$MOSTA.1MODEL1!="" & df2$MEAN_FST_100<=0.3]
  #           ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!="" & df2$MEAN_FST_100<=0.3],
  #           ylab="Mean Fst",col=cols)
  #   boxplot(df2$MEAN_FST_75[df2$MOSTA.1MODEL1!="" & df2$MEAN_FST_75<=0.3]
  #           ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!="" & df2$MEAN_FST_75<=0.3],
  #           ylab="",col=cols,main=spp)
  #   boxplot(df2$MEAN_FST_50[df2$MOSTA.1MODEL1!="" & df2$MEAN_FST_50<=0.3]
  #           ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!="" & df2$MEAN_FST_50<=0.3],
  #           ylab="",col=cols)
  #   mod100=aov(df2$MEAN_FST_100[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
  #   mod75=aov(df2$MEAN_FST_75[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
  #   mod50=aov(df2$MEAN_FST_50[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
  #   print(summary(mod100))
  #   print(summary(mod75))
  #   print(summary(mod50))
  #   print(TukeyHSD(mod100))
  #   print(TukeyHSD(mod75))
  #   print(TukeyHSD(mod50))
    
    print("DXY")
    par(mfrow=c(1,3))
    boxplot(df2$MEAN_DXY_100[df2$MOSTA.1MODEL1!=""]
            ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""],
            ylab="Mean DXY",col=cols)
    boxplot(df2$MEAN_DXY_75[df2$MOSTA.1MODEL1!=""]
            ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""],
            ylab="",col=cols,main=spp)
    boxplot(df2$MEAN_DXY_50[df2$MOSTA.1MODEL1!=""]
            ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""],
            ylab="",col=cols)
    mod100=aov(df2$MEAN_DXY_100[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
    mod75=aov(df2$MEAN_DXY_75[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
    mod50=aov(df2$MEAN_DXY_50[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
    print(summary(mod100))
    print(summary(mod75))
    print(summary(mod50))
    print(TukeyHSD(mod100))
    print(TukeyHSD(mod75))
    print(TukeyHSD(mod50))
    
    # print("RECOMB")
    # par(mfrow=c(1,1))
    # boxplot(df2$MEAN_RECOMB_100[df2$MOSTA.1MODEL1!=""]
    #         ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""],
    #         ylab="Mean Recombination",col=cols)
    # mod100=aov(df2$MEAN_RECOMB_100[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
    # print(summary(mod100))
    # print(TukeyHSD(mod100))

    # print("MISSING")
    # par(mfrow=c(1,1))
    # boxplot(df2$PERCENT_MISSING[df2$MOSTA.1MODEL1!=""]
    #         ~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""],
    #         ylab="Mean Percent Missing",col=cols)
    # mod100=aov(df2$PERCENT_MISSING[df2$MOSTA.1MODEL1!=""]~df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
    # print(summary(mod100))
    # print(TukeyHSD(mod100))
    # 
    # tab = table(df2$MOSTA.1MODEL1[df2$MOSTA.1MODEL1!=""])
    # pie(tab,col=cols,main=spp)
    
  #dev.off()
}

