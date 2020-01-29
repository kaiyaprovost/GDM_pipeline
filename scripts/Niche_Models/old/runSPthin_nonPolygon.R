dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) return(TRUE)
  
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}

packages = c("spThin")

for(p in packages){
  dynamic_require(p)
}

wd = "/Users/kprovost/Documents/Dissertation/GBIF/"
setwd(wd)

filename = "AllSpecies_NotThinned_NoUncertains.txt"
species = "Vireo bellii"

time = as.character(Sys.time())
time = gsub(":","_",time)
time = gsub(" ","_",time)
time = gsub("-","_",time)

logname = paste(species,substring(filename,0,nchar(filename)-4),"_log_file",time,".log",sep="")
outbase = paste(species,substring(filename,0,nchar(filename)-4),"_out_",time,sep="")
outdir = paste(outbase,"/",sep="")
print(c(logname,outbase,outdir))

data = read.csv(filename,header=T,sep="\t",row.names=NULL)
data_small = data[data$species==species,]
print(names(data_small))

## ONLY TAKE DATA WHERE UNCERTAINTY IS LESS THAN 2000

thin(loc.data = data_small, 
    lat.col = "decimallatitude", long.col = "decimallongitude", 
    spec.col = "species", 
    thin.par = 10, reps = 1, 
    locs.thinned.list.return = TRUE, 
    write.files = TRUE, 
    max.files = 1, 
    out.dir = outdir, out.base = outbase, 
    write.log.file = TRUE,
    log.file = logname )
