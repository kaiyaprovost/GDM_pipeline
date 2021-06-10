## DROP THE MIX IN THE UNI AND BIV


#outfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_splines.alltogether"
outfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/bru_20may2021.alltogether"

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")

files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/",
                   pattern="modelparameters",recursive=T,full.names = T)
#files = files[!(grepl("chr",files))]
#files = files[!(grepl("gene",files))]
files = files[(grepl("BRU",files))]

newdf=NULL

for (i in 1:length(files)) {
  print(paste(i,"/",length(files)))
  file = files[i]
  base = tools::file_path_sans_ext(basename(file))
  #print(base)
  types=strsplit(base,"_")[[1]]
  spp=types[1]
  variables=types[2]
  dataset=types[4]
  
  txt=readLines(file,n=20)
  null=txt[grepl("NULL Deviance",txt)]
  dev= txt[grepl("GDM Deviance",txt)]
  exp= txt[grepl("Percent Deviance Explained",txt)]
  
  null=as.numeric(strsplit(null,":")[[1]][2])
  dev= as.numeric(strsplit(dev, ":")[[1]][2])
  exp= as.numeric(strsplit(exp, ":")[[1]][2])
  
  linetoadd = c(spp,dataset,variables,null,dev,exp)
  names(linetoadd) = c("SPECIES","DATASET","MODEL","NULL","DEV","EXP")
  
  if(is.null(newdf)){
    newdf = t(as.data.frame(linetoadd))
  } else {
    newdf = rbind(newdf,linetoadd)
  }
}

newdf=unique(newdf)

write.table(newdf,outfile,row.names = F,append = T)

tab=table(newdf[,"SPECIES"],newdf[,"DATASET"])
corrplot::corrplot(tab,is.corr=F,method="number")


## can you get relative impacts of everything too? 
## get spline heights 
files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/",
                   pattern="splineheights.*txt",recursive=T,full.names = T)

fulldf = data.frame()
for(i in 1:length(files)) {
  print(paste(i,"/",length(files)))
  file = files[i]
  df = read.table(file,header=T,sep=",",colClasses = c("character","numeric"),row.names = 1)
  df = t(df)
  df=as.data.frame(df)
  
  filesplit=strsplit(basename(file),"_")[[1]]
  species=filesplit[1]
  params=filesplit[2]
  chrom=tools::file_path_sans_ext(filesplit[4])
  
  df$SPP = species
  df$PARAMS = params
  df$CHROM = chrom
  
  write.table(df,"~/splineheights.txt",append=T)
  
  fulldf=plyr::rbind.fill(fulldf,df)
  
}

okay=c("ABUN","ENV","IBD","LGM","PRES","STR","STR-ABUN","STR-ENV","STR-IBD","STR-LGM","STR-PRES","STR-IBD-ABUN","STR-IBD-ENV","STR-IBD-LGM","STR-IBD-PRES")

fulldf=fulldf[,c("SPP","PARAMS","CHROM","ABUN","ENV","IBD","LGM","PRES","STR")]
fulldf=unique(fulldf)
fulldf=fulldf[rowSums(is.na(fulldf))!=6,]

fulldf=fulldf[!(grepl("PC",fulldf$PARAMS)),]
fulldf=fulldf[fulldf$PARAMS %in% okay,]

write.table(fulldf,"~/splineheights.txt",append=F)

fulldfna = fulldf
fulldfna[is.na(fulldfna)] = 0
fulldfna$sum = rowSums(fulldfna[,c(4:9)])


fulldf_relative = fulldfna

for(row in 1:nrow(fulldfna)) {
  print(row)
  fulldf_relative[row,4:9] = (fulldfna[row,4:9] / fulldfna[row,10])
}


abun = aggregate(fulldf_relative$ABUN~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(abun) = c("PARAMS","ABUN")
env = aggregate(fulldf_relative$ENV~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(env) = c("PARAMS","ENV")
ibd = aggregate(fulldf_relative$IBD~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(ibd) = c("PARAMS","IBD")
lgm = aggregate(fulldf_relative$LGM~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(lgm) = c("PARAMS","LGM")
pres = aggregate(fulldf_relative$PRES~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(pres) = c("PARAMS","PRES")
str = aggregate(fulldf_relative$STR~fulldf_relative$PARAMS,FUN=function(x){mean(x,na.rm=T)}); colnames(str) = c("PARAMS","STR")

forbarplot = merge(merge(abun,env,all=T),merge(merge(ibd,lgm,all=T),merge(pres,str,all=T),all=T),all=T)
barplot(forbarplot)

par(mar=c(8,6,2,2))
boxplot(fulldf$ABUN~fulldf$PARAMS,las=2,xlab="")
