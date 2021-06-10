## remove old objects and set virtual memory size to handle large vcfs
rm(list=ls())
Sys.setenv('R_MAX_VSIZE'=32000000000)

## function to install or load the needed packages
dynamic_require <- function(package,lib=NULL) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  if(is.null(lib)) {
    install.packages(package,repos='http://cran.us.r-project.org')
  } else {
    install.packages(package,lib=lib,repos='http://cran.us.r-project.org')
  }
  return(eval(parse(text = paste(
    "require(", package,  ")"))))
}

packages = c("ape","vcfR","adegenet","poppr",
             "phangorn","plotrix","rgdal","RColorBrewer","tools","R.utils")

for (p in packages) { dynamic_require(p) }

## set path 
mypath="~/GDMs/"
setwd(mypath)
overwrite=F

## get files
filelist=list.files(mypath,pattern=".vcf$",full.names = T,recursive=F)

## sort files by size
x <- file.info(filelist)
filelist = filelist[match(1:length(filelist),rank(x$size))]
filelist = filelist[complete.cases(filelist)]

## iterate through files
for (vcffile in filelist) {
  print(vcffile)
  
  try({
    
    ## read in vcf
    vcf <- vcfR::read.vcfR(vcffile, verbose = TRUE,limit=1e08)
    ## convert to genlight format
    x <- vcfR::vcfR2genlight(vcf) 
    
    ## remove the vcf from R memory
    rm(vcf)
    
    ## if the ploidy is not set, set it to diploid
    if(max(ploidy(x))!=min(ploidy(x))){
      ploidy(x) = 2
    }
    
    ## calculate the genetic distances, convert to matrix, and output
    x.dist3 <- poppr::bitwise.dist(x,missing_match = T,percent=F,scale_missing = T)
    mat3=as.matrix(x.dist3)
    write.csv(mat3,paste(vcffile,"_distancematrix.csv",sep = ""),quote=F,overwrite=overwrite)
    
  })
}
