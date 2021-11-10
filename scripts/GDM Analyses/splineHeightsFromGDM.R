## DROP THE MIX IN THE UNI AND BIV


#outfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_splines.alltogether"
outfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/temp.alltogether2.temp"

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")

files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/",
                   pattern="modelparameters",recursive=T,full.names = T)
files = files[!(grepl("BAD",files))]
files = files[!(grepl("OLD",files))]
files = files[!(grepl("gz",files))]

newdf=NULL

for (i in 1:length(files)) {
  print(paste(i,"/",length(files)))
  file = files[i]
  date = file.info(file)$ctime
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
  
  linetoadd = c(spp,dataset,variables,null,dev,exp,date)
  names(linetoadd) = c("SPECIES","DATASET","MODEL","NULL","DEV","EXP","DATE")
  
  if(is.null(newdf)){
    newdf = t(as.data.frame(linetoadd))
  } else {
    newdf = rbind(newdf,linetoadd)
  }
}

newdf=unique(newdf)

write.table(newdf,outfile,row.names = F,sep="\t")

tab=table(newdf[,"SPECIES"],newdf[,"DATASET"])
corrplot::corrplot(tab,is.corr=F,method="color")
