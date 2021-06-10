## extract information from the model GDM parameter outputs
setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")

listfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/",pattern="modelparameters",recursive = T,full.names = T)
listnames = c("DATASET","SPECIES","NULL","DEVIANCE","EXPLAINED","PARAMETERS")
df = c()

for(filename in listfiles) {
print(filename)
  
#filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BELLII_STR-ABUN_modelparameters_chr1B.txt"
tosplit=strsplit(filename,"/")[[1]]
split = strsplit(tosplit[length(tosplit)],"_")[[1]]
species=split[1]
params=split[2]
chrom=strsplit(split[4],"\\.")[[1]][1]

lines = readLines(filename)
null = as.numeric(strsplit(lines[14],":  ")[[1]][2])
devi = as.numeric(strsplit(lines[15],":  ")[[1]][2])
expl = as.numeric(strsplit(lines[16],":  ")[[1]][2])

temp = cbind(chrom,species,null,devi,expl,params)

if(is.null(df)){
  df = temp
} else {
  df = rbind(df,temp)
}

}

colnames(df) = listnames
write.csv(df,file="extracted_model_parameters_feb2020.csv",row.names = F)
