#filename = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/BELLII/AddedPredThin_BestModel_Vireo bellii_addedLayers_LGM.asc_-5.38855934143066_-61.0209884643555.rescaled.asc.MORPH-AND-GENE-DISTANCES.txt.converted.txt"

filelist = list.files("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED",
                      pattern="converted.txt$",recursive = T,full.names = T)

for (filename in rev(filelist)) {
  
  print(filename)

cat = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/CATALOG_BY_SPECIES.txt"
split = strsplit(filename,"/")[[1]]
species = split[length(split)-1]

catdf = read.csv(cat,sep="\t")
subcat = catdf[catdf$SPP==species,]

filedf = read.csv(filename,header=F,stringsAsFactors = F)
if(ncol(filedf) == 1) {filedf = read.csv(filename,header=F,stringsAsFactors = F,sep="\t")}
colnums = filedf[1,]
rownums = filedf[,1]

colnums = gsub(" ","",colnums)
colnums = gsub("-","",colnums)
colnums = gsub("\\.","",colnums)
colnums = gsub("SKIN","",colnums)
colnums = gsub("SW18","J",colnums)
colnums = gsub("J0","J",colnums)
colnums = gsub("J0","J",colnums)
colnums = gsub("J0","J",colnums)

rownums = gsub(" ","",rownums)
rownums = gsub("-","",rownums)
rownums = gsub("\\.","",rownums)
rownums = gsub("SKIN","",rownums)
rownums = gsub("SW18","J",rownums)
rownums = gsub("J0","J",rownums)
rownums = gsub("J0","J",rownums)
rownums = gsub("J0","J",rownums)

colvals = (lapply(colnums,FUN=function(x){strsplit(x,"/")[[1]]}))
rowvals = (lapply(rownums,FUN=function(x){strsplit(x,"/")[[1]]}))

tokeep=c()
for(i in 1:length(colvals)) {
  test = colvals[[i]]
  x=which(test %in% catdf$CATALOG.NUMBER)
  if(length(x) > 0){
    tokeep = c(tokeep,i)
  }
}
if(!identical(colvals,rowvals)){
  tokeepr=c()
  for(i in 1:length(rowvals)) {
    test = rowvals[[i]]
    x=which(test %in% catdf$CATALOG.NUMBER)
    if(length(x) > 0){
      tokeepr = c(tokeepr,i)
    }
  }
} else{tokeepr=tokeep}

tokeepr = c(1,tokeepr)
tokeep = c(1,tokeep)
subsetdf = filedf[tokeepr,tokeep]
colnames(subsetdf) = subsetdf[1,]
subsetdf = subsetdf[2:nrow(subsetdf),]

write.csv(subsetdf,paste(filename,"REDUCED.csv",sep=""),row.names = F)

}

