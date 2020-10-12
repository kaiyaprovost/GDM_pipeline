## DROP THE MIX IN THE UNI AND BIV


#outfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/multivariate/extracting_splines.alltogether"
outfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/lostruct_gdms/lostruct_splines_2.alltogether"

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/lostruct_gdms/")

files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/lostruct_gdms/",
                   pattern="modelparameters",recursive=T,full.names = T)


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

write.table(newdf,outfile,row.names = F)

tab=table(newdf[,"SPECIES"],newdf[,"DATASET"])
corrplot::corrplot(tab,is.corr=F,method="color")
