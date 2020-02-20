outfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_best_column.alltogether"

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")

files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate",
                   pattern="variable_deviance",recursive=T,full.names = T)
## test is morph only

#files = files[which(grepl("MORPH",files))]
#files = files[which(!grepl("PC",files))]
#files = files[which(grepl("gene.csv",files))]
#files = files[which(!grepl("SINUA",files))]
#files = files[which(!grepl("1B",files))]


newdf = c()

for (file in files) {
  print(file)
  
  df = read.csv(file,row.names = 1)
  
  if(sum(colSums(!is.na(df))) != 0) {
    
    rightcol=max(which(df[1,] == min(df[1,],na.rm=T)))
    values=(df[,names(df)[rightcol]])
    names(values) = (rownames(df))
    
    if(rightcol==13) {
      toadd = "_UNI"
      rightcol=12
    } else {
      toadd = ""
    }
    
    #rownames(otherdf1)
    rnam=c("PC1T","PC2T","PC3T","PC1P","PC2P","PC3P",
               "IBD","LGM","PRES","ENV","MID","STR","ABUN")
    
    otherdf1 = read.csv(sub("deviance","importance",file)) ## before all three had row.names=1
    otherdf2 = read.csv(sub("deviance","significance",file))
    otherdf3 = read.csv(sub("deviance","permutations",file))
    
    if(nrow(otherdf1) == length(rnam)) {
      rownames(otherdf1) = rnam
      rownames(otherdf2) = rnam
      rownames(otherdf2) = rnam
    } else {
      rownames(otherdf1) = make.unique(as.character(rownames(otherdf1)))
      rownames(otherdf2) = make.unique(as.character(rownames(otherdf2)))
      rownames(otherdf2) = make.unique(as.character(rownames(otherdf3)))
    }
    
    otherdf1=otherdf1[,-1]
    otherdf2=otherdf2[,-1]
    otherdf3=otherdf3[,-1]
    
    val1 = cbind(otherdf1[,names(otherdf1)[min(rightcol,ncol(otherdf1),na.rm=T)]])
    names(val1) = paste(rownames(otherdf1),"imp")
    val2 = cbind(otherdf1[,names(otherdf2)[min(rightcol,ncol(otherdf2),na.rm=T)]])
    names(val2) = paste(rownames(otherdf2),"sig")
    val3 = cbind(otherdf3[,names(otherdf3)[min(rightcol,ncol(otherdf3),na.rm=T)]])
    names(val3) = paste(rownames(otherdf3),"per")
    
    toname = strsplit(file,"/")[[1]]
    
    if (is.null(newdf)) {
      
      
      
      newdf = t(as.data.frame(c(paste(toname[length(toname)],toadd,sep=""),
                                values,val1,val2,val3)))
    } else {
      newrow = t(as.data.frame(c(paste(toname[length(toname)],toadd,sep=""),
                                 values,val1,val2,val3)))
      newdf = rbind(newdf,newrow)
    }
    
    rownames(newdf) = NULL
    
  }
}
write.table(newdf,outfile,append=F,row.names = F,sep="\t")




newdf=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_best_column.alltogether",
                 sep="\t",header=T)

## imp cols = 6:18, sig cols == 19:30, per cols = 31:44
## metadata 45:48
imps = as.data.frame(newdf[,c(1:5,6:18,45:48)])
sigs = as.data.frame(newdf[,c(1:5,19:30,45:48)])
#tospp = as.character(imps[,1])
#imps$SPECIES = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][1]}))
#sigs$SPECIES = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][1]}))

#imps$VARIABLE = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][4]}))
#sigs$VARIABLE = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][4]}))
imps$VARIABLE = imps$CHROM
sigs$VARIABLE = sigs$CHROM

imps = imps[order(imps$VARIABLE),]

impsnums = imps[,8:18]
impsnums = as.matrix(impsnums)


#impsnums = impsnums[1:25,]

for(i in 1:nrow(impsnums)) {
  for (j in 1:ncol(impsnums)) {
    
    if(is.na(impsnums[i,j])) {
      #print("NA")
      impsnums[i,j] = 0
    }
    
    if(is.infinite(impsnums[i,j])) {
      #print("inf")
      impsnums[i,j] = 0
    }
    
    if(impsnums[i,j] < 0) {
      #print("NEG")
      impsnums[i,j] = 0
    }
    
    impsnums[i,j] = as.numeric(impsnums[i,j])
  }
  
}

#barplot(t(impsnums),ylim=c(0,100))

## this gets shitty 
for(i in 1:nrow(impsnums)) {
  explain = as.numeric(as.character(imps$`Percent deviance explained`[i]))
  impsnums[i,] = scales::rescale(as.numeric(impsnums[i,]))
  impsnums[i,] = ((as.numeric(impsnums[i,]) / sum(as.numeric(impsnums[i,]))))
  #print(impsnums[i,])
}

#barplot(t(impsnums))


## get relative explanatory power
for(i in 1:nrow(impsnums)) {
  impsnums[i,] = as.numeric(as.character(impsnums[i,]))*as.numeric(as.character(imps$`Percent deviance explained`[i]))
}

#barplot(t(impsnums))

colnames(impsnums) = colnames(imps[,8:18])
rownames(impsnums) = imps$SPECIES




impsnums=data.matrix(impsnums)
rownames(impsnums)=NULL
colnames(impsnums)=NULL

dimmy=dim(impsnums) #40,11 or 116,11
impsnums=as.numeric(impsnums)
impsnums = matrix(impsnums,ncol=11)

ifremove = rowSums(impsnums,na.rm = T)
ifremove[is.na(ifremove)] = 0
ifremove[is.infinite(ifremove)] = 0
toremove = which(ifremove==0)

if(length(toremove)!= 0) {
  impsnums = impsnums[-toremove,]
  ## remove if entire row is 0
  
  env = rowSums(impsnums[,c(1,2,3,4,8)])
  his = rowSums(impsnums[,c(6,7,9,10)])
  ibd = rowSums(impsnums[,c(5,11)])
  
  combo = cbind(env,his,ibd)
  combo=as.data.frame(combo)
  combo$SPECIES = imps$SPECIES[-toremove]
  combo$VARIABLE = imps$VARIABLE[-toremove]
  
} else {

env = rowSums(impsnums[,c(1,2,3,4,8)])
his = rowSums(impsnums[,c(6,7,9,10)])
ibd = rowSums(impsnums[,c(5,11)])

combo = cbind(env,his,ibd)
combo=as.data.frame(combo)
combo$SPECIES = imps$SPECIES
combo$VARIABLE = imps$VARIABLE

}



short = substr(combo$SPECIES,1,3)
shorta = gsub(".csv","",combo$VARIABLE)
shorta[shorta=="test"] = "allmorph"

#short2 = substr(combo$VARIABLE,4,6)
#short2[short2=="tes"] = "all"
ids = paste(shorta,short,sep=".")

for(spp in unique(combo$SPECIES)) {
  print(spp)
  subcombo = combo[combo$SPECIES==spp,]
  plotmax = round((max(rowSums(subcombo[,1:3]))),-1)
  
  png(paste(spp,"_all_gdm_results_bk.png",sep=""),width=900)
  barplot(t(subcombo[,1:3]),col=c("#d7191c","#ffffbf","#2b83ba"),
          names.arg = subcombo$VARIABLE,ylab="% Variation Explained",
          las=2,cex.names = 1,main=spp,
          ylim=c(0,plotmax+10))
  legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
         fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
         bty="n",box.col="white",border="white")
  dev.off()
  
  png(paste(spp,"_all_gdm_results.png",sep=""),width=900)
  par(bg = rgb(0,0,0,0), fg = 'white',
      col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
      font=2,xpd=T)
  barplot(t(subcombo[,1:3]),col=c("#d7191c","#ffffbf","#2b83ba"),
          names.arg = subcombo$VARIABLE,ylab="% Variation Explained",
          las=2,cex.names = 1,main=spp,
          ylim=c(0,plotmax+10))
  legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
         fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
         bty="n",box.col="white",border="white")
  dev.off()
}




png("alltogether_gdm_results.png",width=900)
par(bg = rgb(0,0,0,0), fg = 'white',
    col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
    font=2,xpd=T)
barplot(t(combo[,1:3]),col=c("#d7191c","#ffffbf","#2b83ba"),
        names.arg = ids,ylab="% Variation Explained",las=2,cex.names = 0.5)
legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
       fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
       bty="n",box.col="white",border="white")
dev.off()


png("morphology_gdm_results_bk.png",width=900)
#par(bg = rgb(0,0,0,0), fg = 'white',
#    col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
#    font=2,xpd=T)
barplot(t(combo[,1:3]),
        col=c("#d7191c","#ffffbf","#2b83ba"),
        names.arg = combo$SPECIES,
        las=2,cex.names=0.5,
        ylab="% Variation Explained"#,
        #ylim=c(0,50)#,
        #space=c(0,0,0,0,0,0,0,0,0,0,1,
        #        0,0,0,0,0,0,0,0,0,1,
        #        0,0,0,0,0,0,0,0,0,1,
        #        0,0,0,0,0,0,0,0,0)
        )
legend("top",legend=c("ENV","HIS","IBD"),ncol=3,
       fill=c("#d7191c","#ffffbf","#2b83ba"),col=c("#d7191c","#ffffbf","#2b83ba"),
       bty="n",box.col="white",border="white")
dev.off()





par(bg = rgb(0,0,0,0), fg = 'white',
    col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
    font=2,xpd=T)
barplot(t(impsnums),las=2,col=as.numeric(as.factor(colnames(imps)[8:18])))



par(bg = rgb(0,0,0,0), fg = 'white',
    col="white",col.axis="white",col.lab="white",col.main="white",col.sub="white",
    font=2,xpd=T)
barplot(as.numeric(as.character(imps$`Percent deviance explained`)),
        col=as.factor(imps$SPECIES)#,
        #density=35,
        #angle=c(0,45,90,-45)[(as.numeric(as.factor(imps$VARIABLE)))]
        )


## with the new data univariates
uni = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/univariate/univariate_gdm_results.csv",
               sep="\t")
uni_trimmed = uni[,c(1:9,16:19,26:32)]
#uni_trimmed[uni_trimmed$Match.!="N/A",]
uni_trimmed$MAX.NOSTR[uni_trimmed$MAX.NOSTR==-99] = NA
uni_trimmed$MIN.NOSTR[uni_trimmed$MIN.NOSTR==-99] = NA

uni_trimmed=uni_trimmed[!(is.na(uni_trimmed$BESTA.NOSTR)),]
uni_trimmed=uni_trimmed[!(is.na(uni_trimmed$MOSTA.NOSTR)),]
uni_trimmed=uni_trimmed[!(is.na(uni_trimmed$MOSTA.1MODEL)),]
uni_trimmed=uni_trimmed[!(is.na(uni_trimmed$LEASTA.NOSTR)),]


newlevels= sort(unique(c(as.character(uni_trimmed$MOSTA.NOSTR),
                         as.character(uni_trimmed$LEASTA.NOSTR),
                         as.character(uni_trimmed$BESTA.NOSTR),
                "ZERO")))

uni_trimmed$MOSTA.1MODEL[(uni_trimmed$MOSTA.1MODEL=="")] = "MIXED"


uni_trimmed$MOSTA.NOSTR <- factor(uni_trimmed$MOSTA.NOSTR, levels=newlevels)
uni_trimmed$LEASTA.NOSTR <- factor(uni_trimmed$LEASTA.NOSTR, levels=newlevels)
uni_trimmed$BESTA.NOSTR <- factor(uni_trimmed$BESTA.NOSTR, levels=newlevels)
uni_trimmed$MOSTA.1MODEL <- factor(uni_trimmed$MOSTA.1MODEL,
                                  levels=unique(c(as.character(uni_trimmed$MOSTA.1MODEL))))


# uni_trimmed=uni_trimmed[(uni_trimmed$BESTA!=""),]
# uni_trimmed=uni_trimmed[(uni_trimmed$MOSTA!=""),]
# uni_trimmed=uni_trimmed[(uni_trimmed$LEASTA!=""),]

uni_trimmed$BESTA.NOSTR[(uni_trimmed$BESTA.NOSTR=="")] = "ZERO"
uni_trimmed$MOSTA.NOSTR[(uni_trimmed$MOSTA.NOSTR=="")] = "ZERO"
uni_trimmed$LEASTA.NOSTR[(uni_trimmed$LEASTA.NOSTR=="")] = "ZERO"

newlevels= newlevels[-which(newlevels=="")]
uni_trimmed$MOSTA.NOSTR <- factor(uni_trimmed$MOSTA.NOSTR, levels=newlevels)
uni_trimmed$LEASTA.NOSTR <- factor(uni_trimmed$LEASTA.NOSTR, levels=newlevels)
uni_trimmed$BESTA.NOSTR <- factor(uni_trimmed$BESTA.NOSTR, levels=newlevels)

barplot(uni_trimmed$MAX.NOSTR,ylim=c(0,100),
        col=as.numeric(as.factor(uni_trimmed$MOSTA.NOSTR)))
barplot(uni_trimmed$MIN.NOSTR,
        col=as.numeric(as.factor(uni_trimmed$LEASTA.NOSTR)))

colors_old=c("black","brown","darkorange","goldenrod","red",
         "orange","cyan","blue","purple","magenta",
         "grey","green","yellow","white")

colors_2=c("black","brown","red","cyan","blue",
         "purple","magenta","yellow","white")

colors=c("black","brown","darkorange3","goldenrod","cyan",
         "blue","purple","darkgrey","magenta",
         "grey","green","red",
         "orange","yellow","white")


colors_1mod = c("black","blue","red","yellow","grey")

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/univariate/")
for (spp in unique(uni_trimmed$SPECIES)) {
  
  palette(colors)
  
  df = uni_trimmed[uni_trimmed$SPECIES==spp,]
  
  namesmostleast = sort(unique(c(as.character(df$MOSTA.NOSTR),as.character(df$LEASTA.NOSTR))))
  namesbest = sort(unique(as.character(df$BESTA.NOSTR)))
  namesleast = sort(unique(c(as.character(df$LEASTA.NOSTR))))
  namesmost = sort(unique(c(as.character(df$MOSTA.NOSTR))))
  
  mostleast = sort(unique(c(df$MOSTA.NOSTR,df$LEASTA.NOSTR)))
  least = sort(unique(c(df$LEASTA.NOSTR)))
  most = sort(unique(c(df$MOSTA.NOSTR)))
  best = sort(unique(df$BESTA.NOSTR))
  
  namesmost1mod = sort(unique(c(as.character(df$MOSTA.1MODEL))))
  most1mod = sort(unique(c(df$MOSTA.1MODEL)))
  
  png(paste(spp,"_percent_explained.png",sep=""))
  par(#mfrow=c(2,1),
      xpd=TRUE,mar=c(3,4,3,0))
  x = barplot(df$MAX.NOSTR,ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,
          col=as.numeric((df$MOSTA.NOSTR)))
  #text(x,df$MAX+0.5,df$MOST,cex=0.5,las=2,srt=90)
  legend("top",legend=namesmost,ncol = 4,cex=0.75,
         bty="n",
         fill=most)
  #y = barplot(df$MIN.NOSTR,ylab=paste(spp,"Deviance"),cex.names=0.75,cex.axis=0.75,
  #        las=2,names=df$DATASET,
  #        col=as.numeric((df$LEASTA.NOSTR)))
  #legend("top",legend=namesleast,ncol = 4,cex=0.75,
  #       bty="n",
  #       fill=least)
  #axis(3,at=y,labels=df$dataset,las=2,cex.axis=0.75)
  #text(y,0,df$LEAST,cex=0.5,srt=90)
  dev.off()
  
  png(paste(spp,"_percent_explained_log.png",sep=""))
  par(mfrow=c(2,1),
      xpd=TRUE,mar=c(3,4,3,0))
  x = barplot(df$MAX.NOSTR,ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,
              col=as.numeric((df$MOSTA.NOSTR)))
  #text(x,df$MAX+0.5,df$MOST,cex=0.5,las=2,srt=90)
  legend("top",legend=namesmost,ncol = 4,cex=0.75,
         bty="n",
         fill=most)
  y = barplot(log(0.000001+df$MIN.NOSTR),ylab=paste(spp,"Log Deviance"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$dataset,
              col=as.numeric((df$LEASTA.NOSTR)))
  legend("top",legend=namesleast,ncol = 4,cex=0.75,
         bty="n",
         fill=least)
  #axis(3,at=y,labels=df$dataset,las=2,cex.axis=0.75)
  #text(y,0,df$LEAST,cex=0.5,srt=90)
  dev.off()
  
  #readline(prompt="Press [enter] to continue")
  png(paste(spp,"_model_fits.png",sep=""))
  plot(df$MAX.NOSTR,df$MIN.NOSTR,pch=22,bg=as.numeric((df$BESTA.NOSTR)),
       ylab=paste(spp,"Deviance"),xlab=paste(spp,"% Explained"))
  legend("topright",legend=namesbest,ncol = 2,cex=0.75,
         fill=best)
  dev.off()
  
  
  
  png(paste(spp,"_percent_explained_log_conflict.png",sep=""),width=480*1.5)
  par(mfrow=c(2,1),xpd=TRUE,mar=c(3,4,3,0))
  x = barplot(df$MAX.NOSTR,ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,
              col=as.numeric((df$MOSTA.NOSTR)))
  barplot(df$MAX.NOSTR,ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
          las=2,names=df$DATASET,density=10,add=T,
          col=as.numeric((df$LEASTA.NOSTR)))
  #text(x,df$MAX+0.5,df$MOST,cex=0.5,las=2,srt=90)
  legend("topleft",legend=namesmost,ncol = 4,cex=0.75,
         bty="n",
         fill=most)
  y = barplot(log(0.000001+df$MIN.NOSTR),ylab=paste(spp,"Log Deviance"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,ylim=c(0,5),
              col=as.numeric((df$LEASTA.NOSTR)))
  barplot(log(0.000001+df$MIN.NOSTR),ylab=paste(spp,"Log Deviance"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,density=10,add=T,ylim=c(0,5),
              col=as.numeric((df$MOSTA.NOSTR)))
  legend("topleft",legend=namesleast,ncol = 4,cex=0.75,
         bty="n",
         fill=least)
  #axis(3,at=y,labels=df$dataset,las=2,cex.axis=0.75)
  #text(y,0,df$LEAST,cex=0.5,srt=90)
  dev.off()
  
  
  palette(colors_1mod)
  
  toplotlog = log10(df$MAX.NOSTR+0.001)
  #toplotlog[toplotlog==-5]=0
  toplotlog=toplotlog+3
  
  png(paste(spp,"_percent_explained_log_1mod.png",sep=""),width=480*1.5)
  par(mfrow=c(1,1),xpd=TRUE,mar=c(5,4,0.5,0))
  x = barplot(toplotlog,ylab=paste(spp,"Log % Explained"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,ylim=c(0,5),yaxt="n",
              col=as.numeric((df$MOSTA.1MODEL)))
  #text(x,df$MAX+0.5,df$MOST,cex=0.5,las=2,srt=90)
  legend("topleft",legend=namesmost1mod,ncol = 4,cex=0.75,
         bty="n",
         fill=most1mod)
  axis(2,at=c(0,1,2,3,4,5),labels=c(0,0.01,0.1,1,10,100),las=2)
  dev.off()
  
  png(paste(spp,"_percent_explained_1mod.png",sep=""),width=480*1.5)
  par(mfrow=c(1,1),xpd=TRUE,mar=c(5,4,0.5,0))
  x = barplot(df$MAX.NOSTR,ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,
              col=as.numeric((df$MOSTA.1MODEL)))
  #text(x,df$MAX+0.5,df$MOST,cex=0.5,las=2,srt=90)
  legend("topleft",legend=namesmost1mod,ncol = 4,cex=0.75,
         bty="n",
         fill=most1mod)
  dev.off()
  
  png(paste(spp,"_percent_explained_1mod_scaled.png",sep=""),width=480*1.5)
  par(mfrow=c(1,1),xpd=TRUE,mar=c(5,4,0.5,0))
  x = barplot(df$MAX.NOSTR,ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
              las=2,names=df$DATASET,ylim=c(0,100),
              col=as.numeric((df$MOSTA.1MODEL)))
  #text(x,df$MAX+0.5,df$MOST,cex=0.5,las=2,srt=90)
  legend("topleft",legend=namesmost1mod,ncol = 4,cex=0.75,
         bty="n",
         fill=most1mod)
  dev.off()
  
}

library(plotrix)
gap.barplot(df$MAX.NOSTR,gap=c(10,70),
            ylab=paste(spp,"% Explained"),cex.names=0.75,cex.axis=0.75,
            las=2,names=df$DATASET,
            col=as.numeric((df$MOSTA.1MODEL)))
       

dx = uni_trimmed[,c(1,2,20,8)]
dx$MAX.NOSTR
dx$SCALEDMAX = scales::rescale(dx$MAX.NOSTR,to=c(100,255))
dtt = table(uni_trimmed[,1:2])
dttn = dtt
dtte = dtt

for(ri in 1:nrow(dtt)) {
  for(ci in 1:ncol(dtt)) {

    variable = rownames(dtt)[ri]
    species = colnames(dtt)[ci]
    
    entry = dx[dx$DATASET==variable,]
    entry = entry[entry$SPECIES==species,]
    entry_1 = as.character(as.factor(entry$MOSTA.1MODEL))
    entry_2 = as.numeric(as.factor(entry$MOSTA.1MODEL))
    entry_3 = as.numeric(as.character(entry$SCALEDMAX))
    if(is.na(entry_3)){
      entry_3 = 0
    }
    dtt[ri,ci] = entry_1
    dttn[ri,ci] = entry_2
    dtte[ri,ci] = entry_3
    
  }
}

#toorderrows=c(36,28,23,16,27,5,9,7,35,25,1,14,31,30,15,8,3,21,20,33,17,32,19,29,10,22,13,11,6,2,4,12,24,26,34,18,37,38,39,40)
toordercols = c(2,3,9,8,6,10,1,7,5,4)

#dtt=dtt[toorderrows,toordercols]
#dttn=dttn[toorderrows,toordercols]
dtt=dtt[,toordercols]
dttn=dttn[,toordercols]
dtte=dtte[,toordercols]

#write.csv(dtt,"dtt.csv")

ff = factor(as.matrix(dtt))
#ff<-factor(as.matrix(dx[2:4]), 
#           levels=c("Done","WIP",""), 
#           labels=c("done","wip","-empty-")
#)
fx<-matrix(as.numeric(ff), ncol=ncol(dtt))
colnames(fx) = colnames(dtt)
rownames(fx) = rownames(dtt)

#fx=fx[,c(2,3,9,8,6,10,1,7,5,4)]

fxd = as.data.frame(fx)

#use labels to assign colors
col<-c("IBA"="black","IBE"="red","IBD"="blue",
       "IBH"="yellow","MIXED"="grey")

brewcols = RColorBrewer::brewer.pal(12,"Dark2")
blue=brewcols[3] 
yellow=brewcols[6]
red=brewcols[4]
purple=brewcols[1]
black="#000000"

imgflip<-function(x) {t(x[nrow(x):1,])}

png("heatmap_gdm_results_all.png",width=900)
par(mfrow=c(1,1),xpd=TRUE,mar=c(5,10,6,1))
image((dttn),col=c(purple,blue,yellow,red,black),
      xaxt="n",yaxt="n")
axis(3,at=seq(0,1,length.out=nrow(dtt)),
     labels=(rownames(fx)),las=2)
axis(2,at=seq(0,1,length.out=ncol(dtt)),
     labels=(colnames(fx)),las=2)
legend(x = 0,y=-0.1,legend=levels(ff),ncol = 5,
       bty="n",
       fill=c(purple,blue,yellow,red,black))
dev.off()

class(dttn) = "matrix"
class(dtte) = "matrix"

fills=c(purple,blue,yellow,red,black)
col2rgb(fills,alpha=T)

fillsf = fills[dttn]
fillsfc = col2rgb(fillsf,alpha=T)

dttef = matrix(dtte)
fillsfc[4,] = floor(dttef[,1])

transparents = rgb(red=fillsfc[1,],green=fillsfc[2,],
    blue=fillsfc[3,],alpha=fillsfc[4,],
    maxColorValue = 255)



# rgba2hex <- function(r,g,b,a) sprintf('#%s',paste(as.hexmode(c(r,g,b,a)),collapse = ''))
# transparents = sapply(1:ncol(fillsfc),FUN = function(x){
#   rgba2hex(fillsfc[1,x],fillsfc[2,x],fillsfc[3,x],fillsfc[4,x])
# })

transp_f = factor(as.matrix(transparents))
imgcolors = as.numeric(transp_f)
fx<-matrix(as.numeric(imgcolors), ncol=ncol(dtt))
colnames(fx) = colnames(dtt)
rownames(fx) = rownames(dtt)

colmat = sapply(c(purple,blue,yellow,red,black),FUN=function(ramp){
  trans = paste(ramp,"64",sep="") ## 32 for 50, 64 for 100
  solid = paste(ramp,"FF",sep="")
  colfunc <- colorRampPalette(c(solid, trans),alpha=T)
  (matrix(colfunc(75), ncol=1))
})

## for each column make a raster and then plot that 
library(grid)

colmatf = factor(colmat)
colcol = as.numeric(colmatf)
colff<-matrix(as.numeric(colcol), ncol=ncol(colmat))
colnames(colff) = colnames(colmat)
rownames(colff) = rownames(colmat)

png("heatmap_gdm_results_all_transparent.png",width=900)
#layout(matrix(1:2,ncol=2), width = c(5,1),height = c(1,1))
layout(t(matrix(c(1,3,1,3,1,2,1,2,1,3,1,3),nrow=2)), width = c(6,1),height = c(1,1))
par(#mfrow=c(1,1),
    xpd=TRUE,mar=c(1,3,6,1))
image(fx,col=levels(transp_f),
      xaxt="n",yaxt="n")
axis(3,at=seq(0,1,length.out=nrow(dtt)),
     labels=(rownames(fx)),las=2)
axis(2,at=seq(0,1,length.out=ncol(dtt)),
     labels=substr(colnames(fx),1,3),las=2)
legend(x = 0,y=-0.1,legend=levels(ff),ncol = 5,
       bty="n",
       fill=c(purple,blue,yellow,red,black))
image(colff,col=levels(colmatf),
      xaxt="n",yaxt="n",main="% Explained")
axis(3,at=seq(0,1,length.out=6),
     labels=round(seq(75,0,length.out=6)),las=2)
axis(2,at=seq(0,1,length.out=5),
     labels=levels(ff),las=2)
dev.off()
##64-FF




fillsx = matrix(fills[dttn],ncol=ncol(dttn),nrow=nrow(dttn))
#50-255

dttnr = raster(t(as.matrix(dttn)))
dtter = raster(t(as.matrix(dtte)))
dttnr

newdf=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_best_column.alltogether",
                 sep="\t",header=T)


dx = newdf[,c(48:50)]
dtt = table(uni_trimmed[,1:2])
#dtt = t(dtt)
dttn = dtt
dtte

for(ri in 1:nrow(dtt)) {
  for(ci in 1:ncol(dtt)) {
    
    variable = rownames(dtt)[ri]
    species = colnames(dtt)[ci]
    
    entry = dx[dx$CHROM==variable,]
    entry = entry[entry$SPECIES==species,]
    entry_1 = as.character(entry$HIGHEST.MODEL)
    entry_2 = as.numeric(entry$HIGHEST.MODEL)
    entry_3 = as.numeric(entry$SCALEDMAX)
    
    if(length(entry_1)!=0){
      dtt[ri,ci] = entry_1
      dttn[ri,ci] = entry_2
      dtte[ri,ci] = entry_3
      if(is.na(entry_3)){
        entry_3 = 0
      }
    } else {
      dtt[ri,ci] = "NONE"
      dttn[ri,ci] = 5
      dtte[ri,ci] = 0
    }

    
  }
}

#toorderrows=c(36,28,23,16,27,5,9,7,35,25,1,14,31,30,15,8,3,21,20,33,17,32,19,29,10,22,13,11,6,2,4,12,24,26,34,18,37,38,39,40)
toordercols = c(2,3,9,8,6,10,1,7,5,4)

dtt=dtt[toorderrows,toordercols]
dttn=dttn[toorderrows,toordercols]


#write.csv(dtt,"dtt.csv")

ff = factor(as.matrix(dtt),
            levels=c("IBD","IBA","IBE","IBH","NONE"))
#ff<-factor(as.matrix(dx[2:4]), 
#           levels=c("Done","WIP",""), 
#           labels=c("done","wip","-empty-")
#)
fx<-matrix(as.numeric(ff), ncol=ncol(dtt))
colnames(fx) = colnames(dtt)
rownames(fx) = rownames(dtt)

#fx=fx[,c(2,3,9,8,6,10,1,7,5,4)]

fxd = as.data.frame(fx)

#use labels to assign colors
col<-c("IBA"="black","IBE"="red","IBD"="blue",
       "IBH"="yellow","MIXED"="grey")

imgflip<-function(x) {t(x[nrow(x):1,])}

png("heatmap_gdm_results_all_multi.png",width=900)
par(mfrow=c(1,1),xpd=TRUE,mar=c(5,10,6,1))
image((dttn),col=c("black","#2b83ba","#d7191c","#ffffbf","white"),
      xaxt="n",yaxt="n")
axis(3,at=seq(0,1,length.out=nrow(dtt)),
     labels=(rownames(fx)),las=2)
axis(2,at=seq(0,1,length.out=ncol(dtt)),
     labels=(colnames(fx)),las=2)
legend(x = 0,y=-0.1,legend=levels(ff),ncol = 5,
       bty="n",
       fill=c("black","#2b83ba","#d7191c","#ffffbf","white"))
dev.off()



image(imgflip(fx),
      
      #breaks=(1:(nlevels(ff)+1))-.5,
      col=col[levels(ff)],
      xaxt="n", yaxt="n"
)
axis(2, at=seq(0,1,length.out=nrow(dx)), labels=rev(paste("Task",dx$Tasks)), las=2)
axis(3, at=seq(0,1,length.out=length(names(dx))-1), labels=names(dx)[-1])     

            

