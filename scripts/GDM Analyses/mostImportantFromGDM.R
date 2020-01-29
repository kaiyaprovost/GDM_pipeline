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


#newdf=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_best_column.morph",
#                 sep="\t")

## imp cols = 6:18, sig cols == 19:30, per cols = 31:44
imps = as.data.frame(newdf[,c(1:5,6:18)])
sigs = as.data.frame(newdf[,c(1:5,19:30)])
tospp = as.character(imps[,1])
imps$SPECIES = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][1]}))
sigs$SPECIES = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][1]}))

imps$VARIABLE = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][4]}))
sigs$VARIABLE = unlist(lapply(tospp,FUN=function(x){strsplit(x,"_")[[1]][4]}))

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
        col=as.factor(imps$SPECIES),
        #density=35,
        #angle=c(0,45,90,-45)[(as.numeric(as.factor(imps$VARIABLE)))]
        )

